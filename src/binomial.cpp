#include <Rcpp.h>
#include <iostream>
#include "binomial.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List binomial_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
               IntegerVector t, IntegerVector R0,
               IntegerVector l_full, IntegerVector r_full, IntegerVector t_full,
               double toler, int max_it) {

    binomial binomial_ob(lambda, l, r, t, R0, l_full, r_full, t_full, toler, max_it);

    binomial_ob.run();

    List out;
    out["llike"] = binomial_ob.llike;
    out["it"] = binomial_ob.it;
    out["lambda"] = binomial_ob.lambda_0;
    out["surv"] = binomial_ob.surv;

    return out;

}

// outer loop
void binomial::run() {
  double old_like = R_NegInf;
  double cond_trans = 0;
  while (it < maxit && !conv) {
      // vector of derivative contributions per participant
      for (int i = 0; i < n_obs; i++) {
          c[i] = surv[right[i]] / surv[left[i]];

          for (int j = left[i]; j < right[i]; j++) {
              cond_trans = lambda_0[j] / (1 - c[i]);
              n_trans[j] += cond_trans;
          }
      }

      // calculate new lambda values
      for (int j = 0; j < n_int - 1; j++) {
          lambda_1[j] = n_trans[j] / risk_0[j];

        if (lambda_0[j] >= 1) {
          lambda_0[j] = 1 - 1e-10;
        }

          surv[j + 1] = surv[j] * (1 - lambda_1[j]);

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
double binomial::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log(surv[left_full[i]] - surv[right_full[i]]) - 
        log(surv[trun_full[i]]);
  }

  return like;
}

