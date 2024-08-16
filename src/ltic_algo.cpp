#include <Rcpp.h>
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
List newton_algorithm(NumericVector lambda, IntegerVector l, IntegerVector r) {

    int n_int = lambda.length();
    int n_obs = l.length();
    std::vector<double> lambda_0 = Rcpp::as< std::vector<double> >(lambda);
    std::vector<int> left = Rcpp::as< std::vector<int> >(l);
    std::vector<int> right = Rcpp::as< std::vector<int> >(r);
    std::vector<double> lambda_1 = lambda_0;
    std::vector<double> cum_lambda(n_int + 1);
    std::vector<double> c(n_obs), deriv(n_obs);
    std::vector<double> deriv_1(n_int), deriv_2(n_int);
    double tol = 1e-8;
    int it = 0;
    bool conv = false;
    double old_lk = R_NegInf;
    double new_lk = 0;

    // initiate cum_lambda
    for (int j = 1; j < n_int + 1; j++){
        cum_lambda[j] = cum_lambda[j - 1] + lambda_0[j - 1];
    }

    // algorithm
    while (it < 1000 && (!conv)) {

        // calculate derivatives
        // vector of derivative contributions per participant
        for (int i = 0; i < n_obs; i++) {
            c[i] = exp(- (cum_lambda[right[i]] - cum_lambda[left[i]]));

            deriv[i] = c[i] / (1 - c[i]);

            // add up contributions from each participant
            for (int j = 0; j < right[i]; j++) {

                if (j < l[i]) {
                    deriv_1[j] -= 1;
                } else if (j >= left[i] && j < right[i]) {
                    deriv_1[j] += deriv[i];
                    deriv_2[j] -= (deriv[i] + deriv[i] * deriv[i]);
                }

            }
        }

        // newton step
        for (int j = 0; j < n_int; j++) {
            lambda_1[j] = lambda_0[j] - deriv_1[j] / deriv_2[j];

            if (lambda_1[j] < 0) {
                lambda_1[j] = 0;
            }

            //conv = conv * (abs(exp(-lambda_0[j]) - exp(-lambda_1[j])) < tol);

            // reset for next iteration
            lambda_0[j] = lambda_1[j];
            cum_lambda[j + 1] = cum_lambda[j] + lambda_0[j];
            deriv_1[j] = 0;
            deriv_2[j] = 0;

        }

        for (int i = 0; i < n_obs; i++) {
            new_lk += log( exp(-cum_lambda[left[i]]) - exp(-cum_lambda[right[i]]));
        }

        conv = new_lk - old_lk < tol;
        old_lk = new_lk;
        new_lk = 0;

        it++;

    }
    
    List res;
    res["lambda_0"] = lambda_0;
    res["lambda_1"] = lambda_1;
    res["it"] = it;
    res["conv"] = conv;
    res["like"] = old_lk;
    return res;
}

// [[Rcpp::export]]
List em_algorithm(NumericVector lambda, IntegerVector l, IntegerVector r, NumericVector R0) {

    int n_int = lambda.length();
    int n_obs = l.length();
    std::vector<int> left = Rcpp::as< std::vector<int> >(l);
    std::vector<int> right = Rcpp::as< std::vector<int> >(r);
    std::vector<double> lambda_0 = Rcpp::as< std::vector <double> >(lambda);
    std::vector<double> lambda_1 = lambda_0;
    std::vector<double> cum_lambda(n_int + 1);
    std::vector<double> c(n_obs), deriv(n_obs);
    std::vector<double> n_trans(n_int), cum_n_trans(n_int + 1), h(n_int);
    std::vector<double> risk_0 = Rcpp::as< std::vector<double> >(R0);
    double cond_trans, tol = 1e-8;
    int it = 0;
    bool conv = false;
    double old_lk = R_NegInf;
    double new_lk = 0;

    // initiate cum_lambda
    for (int j = 1; j < n_int + 1; j++){
        cum_lambda[j] = cum_lambda[j - 1] + lambda_0[j - 1];
    }

    // algorithm loop
    while (it < 1000 && (!conv)) {

        // vector of derivative contributions per participant
        for (int i = 0; i < n_obs; i++) {
            c[i] = exp(- (cum_lambda[right[i]] - cum_lambda[left[i]]));

            for (int j = left[i]; j < right[i]; j++) {
                    cond_trans = (1 - exp(-lambda_0[j])) * exp(- (cum_lambda[j] - cum_lambda[left[i]])) / (1 - c[i]);
                    n_trans[j] += cond_trans;
            }
        }

        // calculate new h values
        for (int j = 0; j < n_int; j++) {
            cum_n_trans[j + 1] = cum_n_trans[j] + n_trans[j];
            h[j] = n_trans[j] / (risk_0[j] - cum_n_trans[j]);

            if (h[j] < 1) {
                lambda_0[j] = - log(1 - h[j]);
            } else {
                lambda_0[j] = 9999;
            }
            cum_lambda[j + 1] = cum_lambda[j] + lambda_0[j];

            //conv = conv * (abs(lambda_0[j] - lambda_1[j]) < tol);

            // reset for next iteration
            n_trans[j] = 0;
            lambda_1[j] = lambda_0[j];

        }

        for (int i = 0; i < n_obs; i++) {
            new_lk += log( exp(-cum_lambda[left[i]]) - exp(-cum_lambda[right[i]]));
        }

        conv = new_lk - old_lk < tol;
        old_lk = new_lk;
        new_lk = 0;

        it++;
    }

    List res;
    res["h"] = h;
    res["lambda_0"] = lambda_0;
    res["lambda_1"] = lambda_1;
    res["it"] = it;
    res["conv"] = conv;
    res["like"] = old_lk;
    return res;

}


