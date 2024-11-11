#include <Rcpp.h>
#include <iostream>
#include "ltic_turnbull.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List ltic_turn_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
            IntegerVector t, double toler, int max_it) {

    ltic_turn ltic_ob(lambda, l, r, t, toler, max_it);

    ltic_ob.run();

    List out;
    out["llike"] = ltic_ob.llike;
    out["it"] = ltic_ob.it;
    out["lambda"] = ltic_ob.cum_lambda;
    out["tol"] = ltic_ob.tol;
    //out["conv"] = ltic_ob.calc_conv();

    return out;

}

// outer loop
void ltic_turn::run() {
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
double ltic_turn::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log( exp(-cum_lambda[left[i]] + cum_lambda[trun[i]]) - 
        exp(-cum_lambda[right[i]] + cum_lambda[trun[i]]));
  }

  return like;
}

// EM algorithm
void ltic_turn::em_algo() {

    int it_em = 0;
    double total_weight = 0;
    while (it_em < 10) {
      // vector of derivative contributions per participant
      calc_weight_sums();

      // calculate weights values
      for (int j = 0; j < n_int; j++) {
        total_weight += s_0[j] * w_sum[j];
      }

      // calculate new s values
      for (int j = 0; j < n_int; j++) {
          s_1[j] = s_0[j] * w_sum[j] / total_weight;
      }
      
      // reset for next iteration
      for (int j = 0; j < n_int; j++) {
          s_0[j] = s_1[j];
          surv[j + 1] = surv[j] - s_1[j];
          w_sum[j] = 0;
      }
      total_weight = 0;
      surv[n_int] = 0.;
      it_em++;
    }
}

void ltic_turn::calc_weight_sums() {

  int l_size, r_size, t_size;

  // first interval
  int size_0 = lr_inv[0].in.size();
  int curr;
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].in[i];
    w_sum[0] = w_sum[0] + (1. / (surv[left[curr]] - surv[right[curr]]));
  }

  for (int i = 0; i < n_obs; i++) {
    w_sum[0] = w_sum[0] + 1. / surv[trun[i]];
  }

  size_0 = lr_inv[0].trun_time.size();
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].trun_time[i];
    w_sum[0] = w_sum[0] - 1. / surv[trun[curr]];
  }

  // loop through the rest
  for (int j = 1; j < n_int; j++) {

    l_size = lr_inv[j].in.size();
    r_size = lr_inv[j].out.size();
    t_size = lr_inv[j].trun_time.size();
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

    // trunc
    for (int i = 0; i < t_size; i++) {
      curr = lr_inv[j].trun_time[i];
      w_sum[j] = w_sum[j] - (1. / surv[trun[curr]]);
    }
  }

}


void ltic_turn::convert_to_haz() {
  for (int j = 1; j < n_int; j++) {
    if (surv[j] >= 0) {
      cum_lambda[j] = -log(surv[j]);
    } else {
      cum_lambda[j] = 9999.;
    }
  }
}

void ltic_turn::convert_to_surv() {
  for (int j = 0; j < n_int; j++) {
    surv[j + 1] = exp(-cum_lambda[j + 1]);
    s_0[j] = surv[j] - surv[j + 1];
  }
}


void ltic_turn::newton_algo() {

    calc_derivs();
    half_steps();

    // reset values
    for (int j = 0; j < n_int; j++) {
        // reset for next iteration
        deriv_1[j] = 0;
        deriv_2[j] = 0;
    }
}

void ltic_turn::calc_derivs() {
    double surv_diff;
    double derv_right, derv_left;

    for (int i = 0; i < n_obs; i++) {

      derv_left = exp(-cum_lambda[left[i]] + cum_lambda[trun[i]]);
      derv_right = exp(-cum_lambda[right[i]] + cum_lambda[trun[i]]);
      surv_diff = derv_left - derv_right;

      deriv_1[trun[i]] += 1.;
      deriv_1[left[i]] += -derv_left / surv_diff;
      deriv_1[right[i]] += derv_right / surv_diff;

      deriv_2[left[i]] += derv_left / surv_diff - derv_left * derv_left / (surv_diff * surv_diff);
      deriv_2[right[i]] += -derv_right / surv_diff - derv_right * derv_right / (surv_diff * surv_diff);
    }
}

// half stepping
void ltic_turn::half_steps() {
    int tries = 0;
    bool inc_lik = false;
    double temp_lk = calc_like();
    double new_lk;
    double alpha = -1;
    
    int n_weight = n_int - 1;
    double w[n_weight];
    double y[n_weight];
    double diff[n_weight];

    for (int j = 0; j < n_weight; j++) {
      if (deriv_2[j + 1] >= 0.) deriv_2[j + 1] = -1e-9;
      //   y[j] = cum_lambda[j + 1];
      //   w[j] = 1e9;
      // } else {
        y[j] = -deriv_1[j + 1] / deriv_2[j + 1] + cum_lambda[j + 1];
        w[j] = deriv_2[j + 1] / 2;
      //}
      //if (y[j] < 0) y[j] = 0;
    }

    /* PAVA algorithm */
    monotoneC(&n_weight, y, w);

    for (int j = 0; j < n_weight; j++) {
      if (y[j] < 0) y[j] = 0.;
    }

    for(int j = 0; j < n_weight; j++){
        diff[j] = y[j] - cum_lambda[j + 1];
    }

    for (int j = 0; j < n_weight; j++) {
      cum_lambda[j + 1] = cum_lambda[j + 1] + diff[j];
    }

    /* Half stepping */
    new_lk = calc_like();
    inc_lik = new_lk >= temp_lk;

    while (tries < 5 && !inc_lik) {
      alpha *= 0.5;

      for (int j = 0; j < n_weight; j++) {
        cum_lambda[j + 1] = cum_lambda[j + 1] + alpha * diff[j];
      }

      new_lk = calc_like();

      tries++;
      inc_lik = new_lk >= temp_lk;
    }
}
