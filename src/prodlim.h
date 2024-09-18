class prodlim{
  public:
    int n_int;
    int n_obs;
    int n_obs_full;

    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> trun;
    std::vector<double> lambda_0;
    std::vector<double> lambda_1;
    std::vector<double> risk_0;
    std::vector<int> left_full;
    std::vector<int> right_full;
    std::vector<int> trun_full;

    std::vector<double> c;
    std::vector<double> cum_lambda;
    std::vector<double> S;
    std::vector<double> p_obs;
    std::vector<double> n_trans, cum_n_trans, h;

    double tol = 1e-8;
    int it = 0;
    bool conv = false;
    bool inc_lik = false;
    double llike = R_NegInf;

    double calc_like();
    void run();


    // constructor
    prodlim(Rcpp::NumericVector lambda, Rcpp::IntegerVector l, 
            Rcpp::IntegerVector r, Rcpp::IntegerVector t, 
            Rcpp::IntegerVector R0, Rcpp::IntegerVector l_full, 
            Rcpp::IntegerVector r_full, Rcpp::IntegerVector t_full) {

      n_int = lambda.length();
      n_obs = l.length();
      n_obs_full = l_full.length();

      left = Rcpp::as< std::vector<int> >(l);
      right = Rcpp::as< std::vector<int> >(r);
      trun = Rcpp::as<std::vector<int> >(t);
      lambda_0 = Rcpp::as< std::vector<double> >(lambda);
      h = Rcpp::as< std::vector<double> >(lambda);
      lambda_1 = lambda_0;
      risk_0 = Rcpp::as< std::vector<double> >(R0);
      left_full = Rcpp::as< std::vector<int> >(l_full);
      right_full = Rcpp::as< std::vector<int> >(r_full);
      trun_full = Rcpp::as<std::vector<int> >(t_full);

      p_obs.resize(n_obs);
      cum_lambda.resize(n_int + 1);
      n_trans.resize(n_int);
      cum_n_trans.resize(n_int + 1);
      //h.resize(n_int, 1/n_int);
      S.resize(n_int + 1, 1);

      // initiate cum_lambda
      for (int j = 1; j < n_int + 1; j++){
          S[j] = S[j - 1] * (1 - h[j - 1]);
      }
      S[n_int] = 0;
      h[n_int - 1] = 1;
      lambda_1[n_int - 1] = R_PosInf;

    };

};
