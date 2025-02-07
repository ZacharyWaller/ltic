#include <Rcpp.h>
#include <iostream>
#include "ltic_s.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List ltic_s_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
            IntegerVector t, IntegerVector R0,  IntegerVector l_full, 
            IntegerVector r_full, IntegerVector t_full, double toler, int max_it) {

    ltic_s ltic_ob(lambda, l, r, t, R0, l_full, r_full, t_full, toler, max_it);

    ltic_ob.run();

    List out;
    out["llike"] = ltic_ob.llike;
    out["it"] = ltic_ob.it;
    out["surv"] = ltic_ob.surv;
    out["tol"] = ltic_ob.tol;

    return out;

}

// outer loop
void ltic_s::run() {
  double old_like = R_NegInf;
  while (it < maxit && !conv) {
    em_algo();
    convert_to_haz();
    newton_algo();
    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
    convert_to_surv();
  }
}

// calculate likelihood
double ltic_s::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log(dist[right_full[i]] - dist[left_full[i]]) - log(1 - dist[trun_full[i]]);
  }

  if (isnan(like)) {
    like = R_NegInf;
  }

  return like;
}

// calculate convergence 
std::vector<double> ltic_s::calc_conv() {
  for (int i = 0; i < n_obs_full; i++) {
    for (int j = trun_full[i]; j < n_obs; j++) {
      deriv_1[j] -= 1 / (surv[trun_full[i]]);
      if (j >= left_full[i] && j < right_full[i]) {
        deriv_1[j] += 1 / (surv[left_full[i]] - surv[right_full[i]]);
      }
    }
  }
  return deriv_1;
}


// EM algorithm
void ltic_s::em_algo() {

    int it_em = 0;
    double cond_trans = 0;
    while (it_em < 10) {
      // vector of derivative contributions per participant
      calc_weight_sums();

      // calculate new h values
      for (int j = 0; j < n_int - 1; j++) {
        cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
        h[j] = n_trans[j] / (risk_0[j] - cum_n_trans[j]);

        if (h[j] >= 1.) {
          h[j] = 1. - 1e-10;
        } else if (isnan(h[j])) {
          h[j] = 0.;
        }

        surv[j + 1] = surv[j] * (1. - h[j]);

        // reset for next iteration
        n_trans[j] = 0.;
        w_sum[j] = 0.;
      }
      it_em++;
    }
}

void ltic_s::calc_weight_sums() {

  int l_size, r_size;

  // first interval
  int size_0 = lr_inv[0].in.size();
  int curr;
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].in[i];
    w_sum[0] = w_sum[0] + (1. / (surv[left[curr]] - surv[right[curr]]));
  }


  // loop through the rest
  for (int j = 1; j < n_int; j++) {

    l_size = lr_inv[j].in.size();
    r_size = lr_inv[j].out.size();
    w_sum[j] = w_sum[j - 1];

    // in
    for (int i = 0; i < l_size; i++) {
      curr = lr_inv[j].in[i];
      w_sum[j] = w_sum[j] + (1. / (surv[left[curr]] - surv[right[curr]]));
    }

    // out
    for (int i = 0; i < r_size; i++) {
      curr = lr_inv[j].out[i];
      w_sum[j] = w_sum[j] - (1. / (surv[left[curr]] - surv[right[curr]]));
    }
  }

  for (int j = 0; j < n_int; j++){
    n_trans[j] = surv[j] * h[j] * w_sum[j];
  }

}

void ltic_s::convert_to_haz() {
  for (int j = 1; j < n_int; j++) {
    dist[j] = 1. - surv[j];
  }
}

void ltic_s::convert_to_surv() {
  for (int j = 0; j < n_int; j++) {
    surv[j + 1] = 1. - dist[j + 1];
    h[j] = 1. - surv[j + 1] / surv[j];
  }
}


void ltic_s::newton_algo() {

    calc_derivs();
    lm_steps();

    // reset values
    for (int j = 0; j < n_int; j++) {
        // reset for next iteration
        deriv_1[j] = 0;
        deriv_2[j] = 0;
    }
}

void ltic_s::calc_derivs() {
    double surv_diff;

    for (int i = 0; i < n_obs_full; i++) {
      surv_diff = surv[left_full[i]] - surv[right_full[i]];

      deriv_1[trun_full[i]] += 1. / surv[trun_full[i]];
      deriv_1[left_full[i]] += -1. / surv_diff;
      deriv_1[right_full[i]] += 1. / surv_diff;

      deriv_2[left_full[i]] += -1. / (surv_diff * surv_diff);
      deriv_2[right_full[i]] += -1. / (surv_diff * surv_diff);
      deriv_2[trun_full[i]] += 1. / (surv[trun_full[i]] * surv[trun_full[i]]);
    }
}

void ltic_s::lm_steps() {

  double temp_llike = R_NegInf;
  int n_weight = n_int - 1;
  int tries = 0;

  for (int j = 0; j < n_weight; j++) {
    temp_dist[j + 1] = dist[j + 1];
  }

  while ((temp_llike < llike || temp_llike == R_NegInf) && tries < 10) {

    ltic_s::icm_step();
    temp_llike = calc_like();

    // reset dist
    if (temp_llike < llike || temp_llike == R_NegInf) {
      for (int j = 0; j < n_weight; j++) {
        dist[j + 1] = temp_dist[j + 1];
      }
      lm_factor *= 0.5;
    } else {
      lm_factor *= 2.;
    }
    tries ++;
  }
}


/* ICM step */
void ltic_s::icm_step() {
    int tries = 0;
    bool inc_lik = false;
    double temp_lk = calc_like();
    double new_lk;
    alpha = -1;
    
    int n_weight = n_int - 1;
    double w[n_weight];
    double y[n_weight];

    /* Newton-like step */
    for (int j = 0; j < n_weight; j++) {
      if (deriv_2[j + 1] >= 0.) deriv_2[j + 1] -= lm_factor;
      y[j] = -deriv_1[j + 1] / deriv_2[j + 1] + dist[j + 1];
      w[j] = deriv_2[j + 1] / 2.;
    }

    /* PAVA algorithm */
    monotoneC(&n_weight, y, w);

    /* Reset any values below 0 or above 1 */
    for (int j = 0; j < n_weight; j++) {
      if (y[j] < 0) y[j] = 0.;
      if (y[j] > 1) y[j] = 1.;
    }

    /* take differences for half-stepping */
    for(int j = 0; j < n_weight; j++){
        diff[j] = y[j] - dist[j + 1];
    }

    for (int j = 0; j < n_weight; j++) {
      dist[j + 1] = dist[j + 1] + diff[j];
    }

    /* Half stepping */
    // new_lk = calc_like();
    // inc_lik = new_lk >= temp_lk;

    // while (tries < 5 && !inc_lik) {
    //   alpha *= 0.5;

    //   for (int j = 0; j < n_weight; j++) {
    //     dist[j + 1] = dist[j + 1] + alpha * diff[j];
    //   }

    //   new_lk = calc_like();

    //   if (isnan(new_lk)) break;

    //   tries++;
    //   inc_lik = new_lk >= temp_lk;
    // }
}
