#include <Rcpp.h>
#include <iostream>
#include "breslow.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List breslow_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
               IntegerVector t, IntegerVector R0,
               IntegerVector l_full, IntegerVector r_full, IntegerVector t_full,
               double toler, int max_it) {

    breslow breslow_ob(lambda, l, r, t, R0, l_full, r_full, t_full, toler, max_it);

    breslow_ob.run();

    List out;
    out["llike"] = breslow_ob.llike;
    out["it"] = breslow_ob.it;
    out["lambda"] = breslow_ob.cum_lambda;

    return out;

}

// outer loop
void breslow::run() {
  double old_like = R_NegInf;
  double cond_trans = 0;
  while (it < maxit && !conv) {
      // vector of derivative contributions per participant
    calc_weight_sums();

      // calculate new lambda values
      for (int j = 0; j < n_int - 1; j++) {
          lambda_1[j] = lambda_0[j] * w_sum[j] / risk_0[j];

          cum_lambda[j + 1] = cum_lambda[j] + lambda_1[j];

          // reset for next iteration
          w_sum[j] = 0.;
          lambda_0[j] = lambda_1[j];
      }
      llike = calc_like();
      conv = llike - old_like < tol && llike - old_like > -tol;
      old_like = llike;
      it++;
  }
}

void breslow::calc_weight_sums() {

  int l_size, r_size;

  // first interval
  int size_0 = lr_inv[0].in.size();
  int curr;
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].in[i];
    w_sum[0] = w_sum[0] + (1. / (1. - exp(-cum_lambda[right[curr]] + cum_lambda[left[curr]])));
  }

  // loop through the rest
  for (int j = 1; j < n_int; j++) {

    l_size = lr_inv[j].in.size();
    r_size = lr_inv[j].out.size();
    w_sum[j] = w_sum[j - 1];

    // in
    for (int i = 0; i < l_size; i++) {
      curr = lr_inv[j].in[i];
      w_sum[j] = w_sum[j] + (1. / (1. - exp(-cum_lambda[right[curr]] + cum_lambda[left[curr]])));
    }

    // out
    for (int i = 0; i < r_size; i++) {
      curr = lr_inv[j].out[i];
      w_sum[j] = w_sum[j] - (1. / (1. - exp(-cum_lambda[right[curr]] + cum_lambda[left[curr]])));
    }
  }

}


// calculate likelihood
double breslow::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log( exp(-cum_lambda[left_full[i]] + cum_lambda[trun_full[i]]) - 
        exp(-cum_lambda[right_full[i]] + cum_lambda[trun_full[i]]));
  }

  return like;
}

