#include <Rcpp.h>
#include <iostream>
#include "icm.h"
#include "monotone.h"
using namespace Rcpp;


// [[Rcpp::export]]
List icm_r(NumericVector lambda, IntegerVector l, IntegerVector r, 
            IntegerVector t, double toler, int max_it) {

    icm ltic_ob(lambda, l, r, t, toler, max_it);

    ltic_ob.run();

    List out;
    out["llike"] = ltic_ob.llike;
    out["it"] = ltic_ob.it;
    out["lambda"] = ltic_ob.cum_lambda;
    out["tol"] = ltic_ob.tol;
    //out["conv"] = ltic_ob.calc_conv();

    return out;

}

// outer loop
void icm::run() {
  double old_like = R_NegInf;
  while (it < maxit && !conv) {
    newton_algo();
    llike = calc_like();
    conv = llike - old_like < tol && llike - old_like > -tol;
    old_like = llike;
    it++;
  }
}

// calculate likelihood
double icm::calc_like() {
  double like = 0;
  for (int i = 0; i < n_obs; i++) {
      like += log( exp(-cum_lambda[left[i]] + cum_lambda[trun[i]]) - 
        exp(-cum_lambda[right[i]] + cum_lambda[trun[i]]));
  }

  return like;
}


void icm::newton_algo() {

    calc_derivs();
    half_steps();

    // reset values
    for (int j = 0; j < n_int; j++) {
        // reset for next iteration
        deriv_1[j] = 0;
        deriv_2[j] = 0;
    }
}

void icm::calc_derivs() {
    double surv_diff;
    double derv_right, derv_left;

    for (int i = 0; i < n_obs; i++) {

      derv_left = exp(-cum_lambda[left[i]] + cum_lambda[trun[i]]);
      derv_right = exp(-cum_lambda[right[i]] + cum_lambda[trun[i]]);
      surv_diff = derv_left - derv_right;

      deriv_1[trun[i]] += 1.;
      deriv_1[left[i]] += -derv_left / surv_diff;
      deriv_1[right[i]] += derv_right / surv_diff;

      deriv_2[left[i]] += derv_left / surv_diff - derv_left * derv_left / (surv_diff * surv_diff);
      deriv_2[right[i]] += -derv_right / surv_diff - derv_right * derv_right / (surv_diff * surv_diff);
    }
}

// half stepping
void icm::half_steps() {
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
      if (deriv_2[j + 1] >= 0.) deriv_2[j + 1] = -1e-9;
      //   y[j] = cum_lambda[j + 1];
      //   w[j] = 1e9;
      // } else {
        y[j] = -deriv_1[j + 1] / deriv_2[j + 1] + cum_lambda[j + 1];
        w[j] = deriv_2[j + 1] / 2;
      //}
      //if (y[j] < 0) y[j] = 0;
    }

    /* PAVA algorithm */
    monotoneC(&n_weight, y, w);

    for (int j = 0; j < n_weight; j++) {
      if (y[j] < 0) y[j] = 0.;
    }

    for(int j = 0; j < n_weight; j++){
        diff[j] = y[j] - cum_lambda[j + 1];
    }

    for (int j = 0; j < n_weight; j++) {
      cum_lambda[j + 1] = cum_lambda[j + 1] + diff[j];
    }

    /* Half stepping */
    new_lk = calc_like();
    inc_lik = new_lk >= temp_lk;

    while (tries < 15 && !inc_lik) {
      alpha *= 0.5;

      for (int j = 0; j < n_weight; j++) {
        cum_lambda[j + 1] = cum_lambda[j + 1] + alpha * diff[j];
      }

      new_lk = calc_like();

      tries++;
      inc_lik = new_lk >= temp_lk;
    }
}
