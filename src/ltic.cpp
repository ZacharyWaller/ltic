#include <Rcpp.h>
#include <iostream>
#include "ltic.h"
using namespace Rcpp;


// [[Rcpp::export]]
List ltic_r(NumericVector lambda, IntegerVector l, IntegerVector r, IntegerVector R0) {

    ltic ltic_ob(lambda, l, r, R0);

    ltic_ob.run();

    List out;
    out["llike"] = ltic_ob.llike;
    out["it"] = ltic_ob.it;
    out["lambda"] = ltic_ob.lambda_0;

    return out;

}

// outer loop
void ltic::run() {
  double old_like = R_NegInf;
  while (it < 1000 && !conv) {
    //em_algo();
    newton_algo();
    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
  }
}

// calculate likelihood
double ltic::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log( exp(-cum_lambda[left[i]]) - exp(-cum_lambda[right[i]]));
  }

  return like;
}

// EM algorithm
void ltic::em_algo() {

    int it_em = 0;
    double cond_trans = 0;
    while (it_em < 10) {
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
      it_em++;
    }
}

void ltic::newton_algo() {

    calc_derivs();
    half_steps();

    // reset values
    for (int j = 0; j < n_int; j++) {
        // reset for next iteration
        deriv_1[j] = 0;
        deriv_2[j] = 0;
        lambda_0[j] = lambda_1[j];
    }
}

void ltic::calc_derivs() {
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
}

// half stepping
void ltic::half_steps() {
    int tries = 0;
    bool inc_lik = false;
    double new_lk = calc_like();
    double alpha = 2;
    while (tries < 3 && !inc_lik) {
      alpha *= 0.5;

      for (int j = 0; j < n_int - 1; j++) {
        lambda_1[j] = lambda_0[j] - alpha * deriv_1[j] / deriv_2[j];

        if (lambda_1[j] < 0) {
          lambda_1[j] = 0;
        }

        cum_lambda[j + 1] = cum_lambda[j] + lambda_1[j];
      }
      cum_lambda[n_int] = R_PosInf;

      new_lk = calc_like();

      tries++;
      inc_lik = new_lk > llike;
    }
}