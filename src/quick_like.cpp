#include <Rcpp.h>

using namespace Rcpp;

// calculate likelihood
// [[Rcpp::export]]
NumericVector calc_like_r(NumericVector step, IntegerVector left, 
                          IntegerVector right, IntegerVector trun, 
                          int n_obs, int n_int) {

  double like = 0;
  std::vector<double> surv(n_int + 1, 1);

  for (int j = 0; j < n_int; j++) {
    surv[j+1] = surv[j] - step[j];
  }
  
  for (int i = 0; i < n_obs; i++) {
    like -= log(surv[left[i]] - surv[right[i]]) - log(surv[trun[i]]);
  }

  NumericVector like_res = wrap(like);
  return like_res;
}

// calculate derivatives
// [[Rcpp::export]]
NumericVector calc_derivs_r(NumericVector step, IntegerVector left, 
                            IntegerVector right, IntegerVector trun, 
                            int n_obs, int n_int) {

  std::vector<double> deriv_1(n_int, 0);
  std::vector<double> surv(n_int + 1, 1);
  double step_tot = 0.;

  for (int j = 0; j < n_int; j++) {
    step_tot += step[j];
  }

  for (int j = 0; j < n_int; j++) {
    surv[j+1] = surv[j] - step[j] / step_tot;
  }

  // calculate derivatives
  for (int i = 0; i < n_obs; i++) {

    for (int j = trun[i]; j < n_int; j++) {
      deriv_1[j] += 1. / surv[trun[i]];
      if (j >= left[i] & j < right[i]) {
        deriv_1[j] -= 1. / (surv[left[i]] - surv[right[i]]);
      }
    }
  }

  NumericVector deriv_r = wrap(deriv_1);
  return deriv_r;
}
