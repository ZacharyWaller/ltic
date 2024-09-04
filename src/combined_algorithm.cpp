#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// Need:
//  - pre-calculate Turnbull intervals and supply to these functions
//  - pre-calculate the indices of the start and end of each observation
// intervals in terms of Turnbull intervals
//  - pre-calculate changes in risk set due to truncation and right-
// censoring

// Algorithm for Newton-like algorithm for left-truncated interval censored data
// TODO pass vectors by reference
// Pre-calculate contributions from right-censoring and left-truncation 
// [[Rcpp::export]]
List combined_algorithm(NumericVector lambda, IntegerVector l, IntegerVector r, 
NumericVector deriv_1_0, IntegerVector R0) {

    int n_int = lambda.length();
    int n_obs = l.length();

    std::vector<int> left = Rcpp::as< std::vector<int> >(l);
    std::vector<int> right = Rcpp::as< std::vector<int> >(r);
    std::vector<double> risk_0 = Rcpp::as< std::vector<double> >(R0);
    std::vector<double> lambda_0 = Rcpp::as< std::vector<double> >(lambda);
    std::vector<double> lambda_1 = lambda_0;

    std::vector<double> c(n_obs), deriv(n_obs);
    std::vector<double> deriv_1 = Rcpp::as< std::vector<double> >(deriv_1_0);
    std::vector<double> deriv_1_init = Rcpp::as< std::vector<double> >(deriv_1_0);
    std::vector<double> deriv_2(n_int);
    std::vector<double> cum_lambda(n_int + 1);
    std::vector<double> n_trans(n_int), cum_n_trans(n_int + 1), h(n_int);
    std::vector<double> ex_lambda_0(n_int);
    std::vector<double> ex_lambda_1(n_int);

    double cond_trans, tol = 1e-8;
    int it = 0;
    bool conv = false;
    bool inc_lik = false;
    double old_lk = R_NegInf;
    double new_lk = 0;
    double alpha = 1;
    int it_em = 0, it_newt = 0, it_big = 0, tries = 0;

    // initiate cum_lambda
    for (int j = 1; j < n_int + 1; j++){
        cum_lambda[j] = cum_lambda[j - 1] + lambda_0[j - 1];
    }
    cum_lambda[n_int] = R_PosInf;
    lambda_0[n_int - 1] = R_PosInf;
    lambda_1[n_int - 1] = R_PosInf;


    // algorithm
    while (it < 1000 && (!conv)) {

      it_em = 0;
      it_newt = 0;
      // EM steps
      while (it_em < 1 && it < 1000 && (!conv)) {

          // vector of derivative contributions per participant
          for (int i = 0; i < n_obs; i++) {
              c[i] = exp(- (cum_lambda[right[i]] - cum_lambda[left[i]]));

              for (int j = left[i]; j < right[i]; j++) {
                      cond_trans = (1 - exp(-lambda_0[j])) * exp(- (cum_lambda[j] - cum_lambda[left[i]])) / (1 - c[i]);
                      n_trans[j] += cond_trans;
              }
          }

          // calculate new h values
          for (int j = 0; j < n_int - 1; j++) {
              cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
              h[j] = n_trans[j] / (risk_0[j] - cum_n_trans[j]);

              if (h[j] < 1) {
                  lambda_0[j] = - log(1 - h[j]);
              } else {
                  lambda_0[j] = 9999;
              }
              cum_lambda[j + 1] = cum_lambda[j] + lambda_0[j];

              // reset for next iteration
              n_trans[j] = 0;
              lambda_1[j] = lambda_0[j];

          }

        for (int i = 0; i < n_obs; i++) {
            new_lk += log( exp(-cum_lambda[left[i]]) - exp(-cum_lambda[right[i]]));
        }

        conv = new_lk - old_lk < tol && new_lk - old_lk > -tol;
        if (new_lk < old_lk) {
          std::cout << "No need to store this string" << std::endl;; 
        }
        old_lk = new_lk;
        new_lk = 0;

          it++;
          it_em++;
      }

      if (conv) {
        std::cout << "converged in EM step" << std::endl;
      }

      while (it_newt < 1 && it < 1000 && (!conv)) {
      // Newton method -------------
      // calculate derivatives
      // vector of derivative contributions per participant
      for (int i = 0; i < n_obs; i++) {
          c[i] = exp(- (cum_lambda[right[i]] - cum_lambda[left[i]]));

          deriv[i] = c[i] / (1 - c[i]);

          // add up contributions from each participant
          for (int j = 0; j < right[i]; j++) {
                if (j < left[i]) {
                    deriv_1[j] -= 1;
                } else {
                    deriv_1[j] += deriv[i];
                    deriv_2[j] += -deriv[i] - deriv[i] * deriv[i];
                }
          }
      }


      // half stepping
      tries = 0;
      new_lk = 0;
      alpha = 2;
      while (tries < 3 && !inc_lik) {
        new_lk = 0;
        alpha *= 0.5;

        for (int j = 0; j < n_int - 1; j++) {
          // newton step
          // if (deriv_2[j] < tol && deriv_2[j] > -tol) {
          //   lambda_1[j] = lambda_0[j];
          // } else if (std::isnan(deriv_2[j])) {
          //   lambda_1[j] = lambda_0[j];
          // } else {
          lambda_1[j] = lambda_0[j] - alpha * deriv_1[j] / deriv_2[j];
          //}

          if (lambda_1[j] < 0) {
            lambda_1[j] = 0;
          }

          cum_lambda[j + 1] = cum_lambda[j] + lambda_1[j];
        }
        cum_lambda[n_int] = R_PosInf;

        for (int i = 0; i < n_obs; i++) {
            new_lk += log( exp(-cum_lambda[left[i]]) - exp(-cum_lambda[right[i]]));
        }

        tries++;
        inc_lik = new_lk > old_lk;
      }

      // reset values
      for (int j = 0; j < n_int; j++) {
          // reset for next iteration
          deriv_1[j] = 0;
          deriv_2[j] = 0;
          lambda_0[j] = lambda_1[j];

      }

      conv = new_lk - old_lk < tol && new_lk - old_lk > -tol;
      old_lk = new_lk;
      new_lk = 0;
      inc_lik = false;
      it++;
      it_newt++;
      }

      it_big++;

    }


    
    List res;
    res["lambda_0"] = lambda_0;
    res["lambda_1"] = lambda_1;
    res["h"] = h;
    res["it"] = it;
    res["it_big"] = it_big;
    res["conv"] = conv;
    res["like"] = old_lk;
    return res;
}
