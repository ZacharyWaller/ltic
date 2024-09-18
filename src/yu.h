class yu{
  public:
    int n_int;
    int n_obs;
    double n;

    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> trun;
    std::vector<double> s_0;
    std::vector<double> s_1;

    std::vector<double> first;
    std::vector<double> second;
    std::vector<double> surv;

    double tol = 1e-8;
    int it = 0;
    bool conv = false;
    bool inc_lik = false;
    double llike = R_NegInf;

    double calc_like();
    void run();


    // constructor
    yu(Rcpp::NumericVector s, Rcpp::IntegerVector l, Rcpp::IntegerVector r, Rcpp::IntegerVector t) {
      n_int = s.length();
      n_obs = l.length();
      n = static_cast<double>(n_obs);

      left = Rcpp::as< std::vector<int> >(l);
      right = Rcpp::as< std::vector<int> >(r);
      trun = Rcpp::as<std::vector<int> >(t);
      s_0 = Rcpp::as< std::vector<double> >(s);
      s_1 = s_0;

      first.resize(n_int);
      second.resize(n_int);
      surv.resize(n_int + 1);

      // initiate surv
      surv[0] = 1;
      for (int j = 0; j < n_int; j++){
          surv[j + 1] = surv[j] - s[j];
      }

    };

};