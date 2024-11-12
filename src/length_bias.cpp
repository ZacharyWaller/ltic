#include <Rcpp.h>
#include <iostream>
#include "length_bias.h"
using namespace Rcpp;


// [[Rcpp::export]]
List length_bias_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
                  IntegerVector t, IntegerVector R0, NumericVector del_t, 
                  double toler, int max_it) {

    length_bias length_bias_ob(lambda, l, r, t, R0, del_t, toler, max_it);

    length_bias_ob.run();

    List out;
    out["llike"] = length_bias_ob.llike;
    out["it"] = length_bias_ob.it;
    out["lambda"] = length_bias_ob.h;

    return out;

}

// outer loop
void length_bias::run() {
  double old_like = R_NegInf;
  double cond_trans;
  while (it < maxit && !conv) {

    // vector of derivative contributions per participant
    calc_weight_sums();

    // calculate new h values
    for (int j = 0; j < n_int - 1; j++) {
        cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
        h[j] = n_trans[j] / (n * (1 - G[j + 1] / G[0]) - cum_n_trans[j]);

        if (h[j] >= 1) {
          h[j] = 1 - 1e-10;
        } else if (isnan(h[j])) {
          h[j] = 0.;
        }
    }

    // reset for next iteration
    for (int j = 0; j < n_int - 1; j++) {
      surv[j + 1] = surv[j] * (1 - h[j]);
      n_trans[j] = 0;
      w_sum[j] = 0;
    }

    for (int j = n_int - 1; j >=0; j--){
        G[j] = G[j + 1] + surv[j] * delta_t[j];
    }

    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
  }
}

void length_bias::calc_weight_sums() {

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
double length_bias::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log(surv[left[i]] - surv[right[i]]) - log(G[0]);
  }

  return like;
}



