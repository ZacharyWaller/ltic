#include <Rcpp.h>
#include <iostream>
#include "turnbull.h"
using namespace Rcpp;


// [[Rcpp::export]]
List turnbull_r(NumericVector s, IntegerVector l, IntegerVector r, IntegerVector t) {

    turnbull turnbull_ob(s, l, r, t);

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
  while (it < 1e5 && !conv) {

      // calculate new M values
      for (int i = 0; i < n_obs; i++) {
          // transitions
          for (int j = left[i]; j < right[i]; j++) {
              M[j] += s_0[j] / (surv[left[i]] - surv[right[i]]);
          }
          // ghosts
          for (int j = 0; j < trun[i]; j++) {
              M[j] += s_0[j] / (surv[trun[i]]);
          }
      }

      // sum up total M
      for (int j = 0; j < n_int; j++) {
          sum_M += M[j];
      }

      // calculate new s values
      for (int j = 0; j < n_int; j++) {
          s_1[j] = M[j] / sum_M;
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
          M[j] = 0;
      }
      sum_M = 0;
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
