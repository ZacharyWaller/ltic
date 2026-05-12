#include <Rcpp.h>
#include <iostream>
#include "ltic.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List ltic_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
            IntegerVector t, IntegerVector R0,  IntegerVector l_full, 
            IntegerVector r_full, IntegerVector t_full, double toler, int max_it) {

    ltic ltic_ob(lambda, l, r, t, R0, l_full, r_full, t_full, toler, max_it);

    ltic_ob.run();

    List out;
    out["llike"] = ltic_ob.llike;
    out["it"] = 11 * ltic_ob.it;
    out["lambda"] = ltic_ob.cum_lambda;
    out["tol"] = ltic_ob.tol;
    out["conv"] = ltic_ob.calc_conv();
    out["deriv"] = ltic_ob.deriv_1;
    // For risk-set calculations
    ltic_ob.calc_weight_sums();
    out["n"] = ltic_ob.n_trans;
    
    return out;

}

// outer loop
void ltic::run() {
  double old_like = R_NegInf;
  while (it < maxit && !conv) {
    em_algo();
    convert_to_haz();
    newton_algo();
    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    if (std::isnan(llike)) {
      break;
    }
    old_like = llike;
    it++;
    convert_to_surv();
  }
}

// calculate likelihood
double ltic::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log( exp(-cum_lambda[left_full[i]] + cum_lambda[trun_full[i]]) - 
        exp(-cum_lambda[right_full[i]] + cum_lambda[trun_full[i]]));
  }
  //if (std::isnan(like)) like = R_NegInf;

  return like;
}

// calculate convergence 
double ltic::calc_conv() {

  calc_derivs();
  double fenchel_1 = 0.;

  for (int j = 0; j < n_int; j++) {
    fenchel_1 += cum_lambda[j] * deriv_1[j];
  }
  return fenchel_1;
}


// EM algorithm
void ltic::em_algo() {

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
        } else if (std::isnan(h[j])) {
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

void ltic::calc_weight_sums() {

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

void ltic::convert_to_haz() {
  for (int j = 1; j < n_int; j++) {
    cum_lambda[j] = -log(surv[j]);
  }
}

void ltic::convert_to_surv() {
  for (int j = 0; j < n_int; j++) {
    surv[j + 1] = exp(-cum_lambda[j + 1]);
    h[j] = 1 - exp(-cum_lambda[j + 1] + cum_lambda[j]);
  }
}


void ltic::newton_algo() {

    calc_derivs();
    half_steps();

    // reset values
    for (int j = 0; j < n_int; j++) {
        // reset for next iteration
        deriv_1[j] = 0;
        deriv_2[j] = 0;
    }
}

void ltic::calc_derivs() {
    double surv_diff;
    double derv_right, derv_left;

    for (int i = 0; i < n_obs_full; i++) {

      derv_left = exp(-cum_lambda[left_full[i]] + cum_lambda[trun_full[i]]);
      derv_right = exp(-cum_lambda[right_full[i]] + cum_lambda[trun_full[i]]);
      surv_diff = derv_left - derv_right;

      deriv_1[trun_full[i]] += 1.;
      deriv_1[left_full[i]] += -derv_left / surv_diff;
      deriv_1[right_full[i]] += derv_right / surv_diff;

      deriv_2[left_full[i]] += derv_left / surv_diff - derv_left * derv_left / (surv_diff * surv_diff);
      deriv_2[right_full[i]] += -derv_right / surv_diff - derv_right * derv_right / (surv_diff * surv_diff);
    }
}

// half stepping
void ltic::half_steps() {
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
      if (deriv_2[j + 1] >= 0.) {
        //Rcpp::Rcout << "LM correction " << j + 1 << std::endl;
        deriv_2[j + 1] = -1e-9;
      }
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
      //Rcpp::Rcout << "Half step!" << std::endl;
      alpha *= 0.5;

      for (int j = 0; j < n_weight; j++) {
        cum_lambda[j + 1] = cum_lambda[j + 1] + alpha * diff[j];
      }

      new_lk = calc_like();

      tries++;
      inc_lik = new_lk >= temp_lk;
    }
}
