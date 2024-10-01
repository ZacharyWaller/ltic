#include <Rcpp.h>
#include <iostream>
#include "yu.h"
using namespace Rcpp;


// [[Rcpp::export]]
List yu_r(NumericVector s, IntegerVector l, IntegerVector r, IntegerVector t,
          double toler, int max_it) {

    yu yu_ob(s, l, r, t, toler, max_it);

    yu_ob.run();

    List out;
    out["llike"] = yu_ob.llike;
    out["it"] = yu_ob.it;
    out["s"] = yu_ob.s_0;

    return out;

}

// outer loop
void yu::run() {
  double old_like = R_NegInf;
  while (it < maxit && !conv) {

    calc_weight_sums();

    // calculate new s values
    for (int j = 0; j < n_int; j++) {
        s_1[j] = s_0[j] * w_sum[j] / n;
        // if (s_1[j] < 0) {
        //   std::cout << "whaet" << std::endl;
        // }
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
        //if (surv[j + 1] < 0) surv[j + 1] = 0.0;
        w_sum[j] = 0;
    }
  }
}

/* Calculate weights of each observation */
void yu::calc_weight_sums() {

  int l_size, r_size, t_size;

  // first interval
  int size_0 = lr_inv[0].in.size();
  int curr;
  w_sum[0] = n_obs;
  for (int i = 0; i < size_0; i++) {
    curr = lr_inv[0].in[i];
    w_sum[0] = w_sum[0] + (1. / (surv[left[curr]] - surv[right[curr]]));
  }

  size_0 = lr_inv[0].trun_time.size();
  for (int i = 0; i < size_0; i++) {
    w_sum[0] = w_sum[0] - 1.;
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
      w_sum[j] = w_sum[j] - 1. / surv[trun[curr]];
    }
  }

}

// calculate likelihood
double yu::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log(surv[left[i]] - surv[right[i]]) - log(surv[trun[i]]);
  }

  return like;
}
