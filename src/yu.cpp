#include <Rcpp.h>
#include <iostream>
#include "yu.h"
using namespace Rcpp;


// [[Rcpp::export]]
List yu_r(NumericVector s, IntegerVector l, IntegerVector r, IntegerVector t) {

    yu yu_ob(s, l, r, t);

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
  while (it < 1e5 && !conv) {

      for (int i = 0; i < n_obs; i++) {
          // transitions
          for (int j = left[i]; j < right[i]; j++) {
              first[j] += 1 / (surv[left[i]] - surv[right[i]]);
          }
          // ghosts
          if (surv[trun[i]] < 1) {
              for (int j = 0; j < n_int; j++) {
                    if (j < trun[i]) {
                        second[j] += 1;
                    } else {
                        second[j] += 1 - 1 / (surv[trun[i]]);
                    }
              }
          }
      }

      // calculate new s values
      for (int j = 0; j < n_int; j++) {
          s_1[j] = s_0[j] * (first[j] + second[j]) / n;
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
          first[j] = 0;
          second[j] = 0;
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
