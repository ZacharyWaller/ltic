#include <Rcpp.h>

using namespace Rcpp;

// Need:
//  - pre-calculate Turnbull intervals and supply to these functions
//  - pre-calculate the indices of the start and end of each observation
// intervals in terms of Turnbull intervals
//  - pre-calculate changes in risk set due to truncation and right-
// censoring

// [[Rcpp::export]]
NumericVector calc_deriv(NumericVector cum_lambda, NumericVector l, NumericVector r, NumericVector q) {

    NumericVector c, deriv;
    NumericVector deriv_1 = 0, deriv_2 = 0;
    int n = l.size();
    /*
    int n_ints = q.size();
    int n_obs = l.size();
    */
   for (int i = 0; i < n; i++) {
        
        c[i] = exp(-cum_lambda[l[i]]);

   }

    /*
    // vector of derivative contributions per participant
    for (int i = 0; i < n_obs; i++) {
        c[i] = exp(-(cum_lambda[l] - cum_lambda[r]));

        deriv[i] = c[i] / (1 - c[i]);

        // add up contributions from each participant
        for (int j = 0; j < n_ints; j++) {

            if (j <= r[i] && j > l[i]) {
                deriv_1[j] += deriv[i];
                deriv_2[j] += - deriv[i] - deriv[i] * deriv[i];
            }

        }
    }
    */
   
   return c;
} 