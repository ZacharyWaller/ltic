#include <Rcpp.h>
#include <iostream>
#include "shen.h"
using namespace Rcpp;


// [[Rcpp::export]]
List shen_r(NumericVector s, IntegerVector l, IntegerVector r, IntegerVector t) {

    shen shen_ob(s, l, r, t);

    shen_ob.run();

    List out;
    out["llike"] = shen_ob.llike;
    out["it"] = shen_ob.it;
    out["s"] = shen_ob.s_0;

    return out;

}

// outer loop
void shen::run() {
  double old_like = R_NegInf;
  while (it < 1e5 && !conv) {

      // calculate new M values
      calc_weight_sums();

      // calculate new s values
      for (int j = 0; j < n_int; j++) {
          s_1[j] = s_0[j] * (1 + w_sum[j] / trun_weight);
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
  }
}

/* Calculate weights of each observation */
void shen::calc_weight_sums() {

  int l_size, r_size, t_size;

  // first interval
  int size_0 = lr_inv[0].in.size();
  int curr;
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].in[i];
    w_sum[0] = w_sum[0] + (1. / (surv[left[curr]] - surv[right[curr]]));
  }

  size_0 = lr_inv[0].trun_time.size();
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].trun_time[i];
    w_sum[0] = w_sum[0] - (1. / surv[trun[curr]]);
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

  trun_weight = 0;
  for (int i = 0; i < n_obs; i++){
    trun_weight += 1. / surv[trun[i]];
  }

}

// calculate likelihood
double shen::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log(surv[left[i]] - surv[right[i]]) - log(surv[trun[i]]);
  }

  return like;
}
