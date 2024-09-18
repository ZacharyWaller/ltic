#include <Rcpp.h>
#include <iostream>
#include "prodlim.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List prodlim_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
               IntegerVector t, IntegerVector R0, IntegerVector l_full,
               IntegerVector r_full, IntegerVector t_full) {

    prodlim prodlim_ob(lambda, l, r, t, R0, l_full, r_full, t_full);

    prodlim_ob.run();

    List out;
    out["llike"] = prodlim_ob.llike;
    out["it"] = prodlim_ob.it;
    out["lambda"] = prodlim_ob.cum_lambda;

    return out;

}

// outer loop
void prodlim::run() {
  double old_like = R_NegInf;
  while (it < 100000 && !conv) {

    double cond_trans = 0;
    // vector of derivative contributions per participant
    for (int i = 0; i < n_obs; i++) {
        p_obs[i] = S[right[i]] - S[left[i]];

        for (int j = left[i]; j < right[i]; j++) {
            cond_trans = h[j] S[j] / p_obs[i];
            n_trans[j] += cond_trans;
        }
    }

    // calculate new h values
    for (int j = 0; j < n_int - 1; j++) {
        cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
        h[j] = n_trans[j] / (risk_0[j] - cum_n_trans[j]);

        if (h[j] < 1) {
            lambda_1[j] = - log(1 - h[j]);
        } else {
            lambda_1[j] = 20;
        }
        cum_lambda[j + 1] = cum_lambda[j] + lambda_1[j];

        // reset for next iteration
        n_trans[j] = 0;
        lambda_0[j] = lambda_1[j];
    }

    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
  }
}

// calculate likelihood
double prodlim::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log(p_obs[left_full[i]] - p_obs[right_full[i]]) - 
        log(S[trun_full[i]]);
  }

  return like;
}



