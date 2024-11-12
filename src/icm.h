class icm{
  public:
    int n_int;
    int n_obs;
    double n;


    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> trun;

    std::vector<double> deriv_1;
    std::vector<double> deriv_2;
    std::vector<double> cum_lambda;

    double tol;
    int maxit;
    int it = 0;
    bool conv = false;
    double llike = R_NegInf;

    double calc_like();
    void newton_algo();
    void calc_derivs();
    void half_steps();
    void run();


    // constructor
    icm(Rcpp::NumericVector lambda, Rcpp::IntegerVector l, 
         Rcpp::IntegerVector r, Rcpp::IntegerVector t,
         double toler, int max_it) {

      tol = toler;
      maxit = max_it;
      n_int = lambda.length();
      n_obs = l.length();
      n = static_cast<double>(n_obs);

      left = Rcpp::as< std::vector<int> >(l);
      right = Rcpp::as< std::vector<int> >(r);
      trun = Rcpp::as<std::vector<int> >(t);

      deriv_1.resize(n_int);
      deriv_2.resize(n_int);
      cum_lambda.resize(n_int + 1);

      /* Initiate Cumulative Hazards */
      for (int j = 1; j < n_int + 1; j++){
          cum_lambda[j] = cum_lambda[j - 1] + lambda[j - 1];
      }
      cum_lambda[0] = 0;
      cum_lambda[n_int] = R_PosInf;

    };

};
