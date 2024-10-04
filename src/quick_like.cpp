#include <Rcpp.h>

using namespace Rcpp;

// calculate likelihood
// [[Rcpp::export]]
NumericVector calc_like_r(NumericVector cum_lambda, IntegerVector left, 
                        IntegerVector right, IntegerVector trun, 
                        int n_obs, int n_int) {

  double like = 0;
  
  for (int i = 0; i < n_obs; i++) {
    like += log( exp(-cum_lambda[left[i]] + cum_lambda[trun[i]]) - exp(-cum_lambda[right[i]] + cum_lambda[trun[i]]));
  }

  NumericVector like_res = wrap(like);
  return like_res;
}

// calculate derivatives
// [[Rcpp::export]]
NumericVector calc_derivs_r(NumericVector cum_lambda, IntegerVector left, 
                                IntegerVector right, IntegerVector trun, 
                                int n_obs, int n_int) {

  double surv_diff;
  double derv_right, derv_left;
  std::vector<double> deriv_1(n_int, 0);

  // calculate derivatives
  for (int i = 0; i < n_obs; i++) {

    derv_left = exp(-cum_lambda[left[i]] + cum_lambda[trun[i]]);
    derv_right = exp(-cum_lambda[right[i]] + cum_lambda[trun[i]]);
    surv_diff = derv_left - derv_right;

    deriv_1[trun[i]] += 1.;
    deriv_1[left[i]] += -derv_left / surv_diff;
    deriv_1[right[i]] += derv_right / surv_diff;
  }

  NumericVector deriv_r = wrap(deriv_1);
  return deriv_r;
}
