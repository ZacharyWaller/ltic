#include <Rcpp.h>
#include <iostream>
#include "turnbull.h"
using namespace Rcpp;


// [[Rcpp::export]]
List turnbull_r(NumericVector s, IntegerVector l, IntegerVector r, IntegerVector t,
                double toler, int max_it) {

    turnbull turnbull_ob(s, l, r, t, toler, max_it);

    turnbull_ob.run();

    List out;
    out["llike"] = turnbull_ob.llike;
    out["it"] = turnbull_ob.it;
    out["s"] = turnbull_ob.s_0;

    return out;

}

// outer loop
void turnbull::run() {
  double old_like = R_NegInf;
  while (it < maxit && !conv) {

      // calculate new M values
      calc_weight_sums();

      for (int j = 0; j < n_int; j++) {
        total_weight += s_0[j] * w_sum[j];
      }

      // calculate new s values
      for (int j = 0; j < n_int; j++) {
          s_1[j] = s_0[j] * w_sum[j] / total_weight;
      }
      
      // check convergence
      llike = calc_like();
      conv = llike - old_like < tol && llike - old_like > -tol;
      old_like = llike;
      it++;

      // reset for next iteration
      for (int j = 0; j < n_int; j++) {
          s_0[j] = s_1[j];
          surv[j + 1] = surv[j] - s_1[j];
          w_sum[j] = 0;
      }
      total_weight = 0;
  }
}

/* Calculate weights of each observation */
void turnbull::calc_weight_sums() {

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

// calculate likelihood
double turnbull::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log(surv[left[i]] - surv[right[i]]) - log(surv[trun[i]]);
  }

  return like;
}
