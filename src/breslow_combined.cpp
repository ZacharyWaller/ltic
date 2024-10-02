#include <Rcpp.h>
#include <iostream>
#include "breslow_combined.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List bres_comb_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
                IntegerVector t, IntegerVector R0, IntegerVector l_full, 
                IntegerVector r_full, IntegerVector t_full, double toler, int max_it) {

    bres_comb b_comb_ob(lambda, l, r, t, R0, l_full, r_full, t_full, toler, max_it);

    b_comb_ob.run();

    List out;
    out["llike"] = b_comb_ob.llike;
    out["it"] = b_comb_ob.it;
    out["lambda"] = b_comb_ob.cum_lambda;
    b_comb_ob.calc_derivs();
    out["derivs"] = b_comb_ob.deriv_1;

    return out;

}

// outer loop
void bres_comb::run() {
  double old_like = R_NegInf;
  while (it < maxit && !conv) {
    em_algo();
    newton_algo();
    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
  }
}

// calculate likelihood
double bres_comb::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs_full; i++) {
      like += log( exp(-cum_lambda[left_full[i]] + cum_lambda[trun_full[i]]) - 
        exp(-cum_lambda[right_full[i]] + cum_lambda[trun_full[i]]));
  }

  return like;
}

// EM algorithm
void bres_comb::em_algo() {

    int it_em = 0;
    double cond_trans = 0;
    while (it_em < 10) {
      // vector of derivative contributions per participant
      for (int i = 0; i < n_obs; i++) {
          c[i] = exp(- (cum_lambda[right[i]] - cum_lambda[left[i]]));

          for (int j = left[i]; j < right[i]; j++) {
              n_trans[j] += lambda_0[j] / (1 - c[i]);
          }
      }

      // calculate new h values
      for (int j = 0; j < n_int - 1; j++) {
          lambda_1[j] = n_trans[j] / risk_0[j];

          cum_lambda[j + 1] = cum_lambda[j] + lambda_1[j];

          // reset for next iteration
          n_trans[j] = 0;
          lambda_0[j] = lambda_1[j];
      }
      it_em++;
    }
}

void bres_comb::newton_algo() {

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

void bres_comb::calc_derivs() {
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
void bres_comb::half_steps() {
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
