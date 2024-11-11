#include <Rcpp.h>
#include <iostream>
#include "prodlim.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List prodlim_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
               IntegerVector t, IntegerVector R0, IntegerVector l_full,
               IntegerVector r_full, IntegerVector t_full, double toler, int max_it) {

    prodlim prodlim_ob(lambda, l, r, t, R0, l_full, r_full, t_full, toler, max_it);

    prodlim_ob.run();

    List out;
    out["llike"] = prodlim_ob.llike;
    out["it"] = prodlim_ob.it;
    out["lambda"] = prodlim_ob.h;

    return out;

}

// outer loop
void prodlim::run() {
  double old_like = R_NegInf;
  double cond_trans;
  while (it < maxit && !conv) {

    // vector of derivative contributions per participant
    calc_weight_sums();

    // calculate new h values
    for (int j = 0; j < n_int - 1; j++) {
        cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
        h[j] = n_trans[j] / (risk_0[j] - cum_n_trans[j]);

        if (h[j] >= 1) {
          h[j] = 1 - 1e-10;
        } else if (isnan(h[j])) {
          h[j] = 0.;
        }

        // reset for next iteration
        surv[j + 1] = surv[j] * (1 - h[j]);
        n_trans[j] = 0;
        w_sum[j] = 0;
    }

    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
  }
}

void prodlim::calc_weight_sums() {

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

// calculate likelihood
double prodlim::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log(surv[left_full[i]] - surv[right_full[i]]) - 
        log(surv[trun_full[i]]);
  }

  return like;
}



