#include <Rcpp.h>
#include <iostream>
#include "ltic.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List ltic_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
            IntegerVector t, IntegerVector R0,  IntegerVector l_full, 
                IntegerVector r_full, IntegerVector t_full) {

    ltic ltic_ob(lambda, l, r, t, R0, l_full, r_full, t_full);

    ltic_ob.run();

    List out;
    out["llike"] = ltic_ob.llike;
    out["it"] = ltic_ob.it;
    out["lambda"] = ltic_ob.cum_lambda;

    return out;

}

// outer loop
void ltic::run() {
  double old_like = R_NegInf;
  while (it < 100000 && !conv) {
    em_algo();
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
  for (int i = 0; i < n_obs_full; i++) {
      like += log( exp(-cum_lambda[left_full[i]] + cum_lambda[trun_full[i]]) - 
        exp(-cum_lambda[right_full[i]] + cum_lambda[trun_full[i]]));
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
              n_trans[j] += (1 - exp(-lambda_0[j])) * exp(- (cum_lambda[j] - cum_lambda[left[i]])) / (1 - c[i]);
          }
      }

      // calculate new h values
      for (int j = 0; j < n_int - 1; j++) {
          cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
          h[j] = n_trans[j] / (risk_0[j] - cum_n_trans[j]);

          if (h[j] < 1) {
              lambda_1[j] = - log(1 - h[j]);
          } else {
              lambda_1[j] = 9999;
          }
          cum_lambda[j + 1] = cum_lambda[j] + lambda_1[j];

          // reset for next iteration
          n_trans[j] = 0;
          lambda_0[j] = lambda_1[j];
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
        lambda_0[j] = cum_lambda[j + 1] - cum_lambda[j];
    }
}

void ltic::calc_derivs() {
    double surv_diff;
    double derv_right, derv_left;

    for (int i = 0; i < n_obs_full; i++) {

      derv_left = exp(-cum_lambda[left_full[i]] + cum_lambda[trun_full[i]]);
      derv_right = exp(-cum_lambda[right_full[i]] + cum_lambda[trun_full[i]]);
      surv_diff = derv_left - derv_right;

      deriv_1[trun_full[i]] += 1;
      deriv_1[left_full[i]] += -derv_left / surv_diff;
      deriv_1[right_full[i]] += derv_right / surv_diff;

      deriv_2[left_full[i]] += derv_left / surv_diff - derv_left * derv_left / (surv_diff * surv_diff);
      deriv_2[right_full[i]] += -derv_right / surv_diff - derv_right * derv_right / (surv_diff * surv_diff);
    }
}

// half stepping
void ltic::half_steps() {
    int tries = 0;
    bool inc_lik = false;
    double temp_lk = calc_like();
    double new_lk;
    double alpha = -1;
    
    int n_weight = n_int - 1;
    double w[n_weight];
    double y[n_weight];
    double diff[n_weight];

    for (int j = 0; j < n_weight; j++) {
      w[j] = deriv_2[j + 1] / 2;
      y[j] = -deriv_1[j + 1] / deriv_2[j + 1] + cum_lambda[j + 1];
    }

    
    monotoneC(&n_weight, y, w);

    for(int j = 0; j < n_weight; j++){
        diff[j] = y[j] - cum_lambda[j + 1];
    }

    for (int j = 0; j < n_weight; j++) {
      cum_lambda[j + 1] = cum_lambda[j + 1] + diff[j];
    }
    new_lk = calc_like();
    inc_lik = new_lk >= temp_lk;

    while (tries < 5 && !inc_lik) {
      alpha *= 0.5;

      for (int j = 0; j < n_weight; j++) {
        cum_lambda[j] = cum_lambda[j] + alpha * diff[j];
      }
      new_lk = calc_like();

      tries++;
      inc_lik = new_lk >= temp_lk;
    }
}
