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
      for (int i = 0; i < n_obs; i++) {
          // transitions
          for (int j = left[i]; j < right[i]; j++) {
              numer[j] += s_0[j] / (surv[left[i]] - surv[right[i]]);
          }
          // ghosts
          for (int j = trun[i]; j < n_int; j++) {
              numer[j] -= s_0[j] / surv[trun[i]];
          }
      }

      // sum up total M
      for (int i = 0; i < n_obs; i++) {
          denom += 1 / surv[trun[i]];;
      }

      // calculate new s values
      for (int j = 0; j < n_int; j++) {
          s_1[j] = s_0[j] + numer[j] / denom;
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
          numer[j] = 0;
      }
      denom = 0;
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
