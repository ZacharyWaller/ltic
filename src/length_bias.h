class length_bias{
  public:
    int n_int;
    int n_obs;
    double n;

    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> trun;
    std::vector<double> risk_0;
    std::vector<double> delta_t;

    std::vector<double> surv;
    std::vector<double> G;
    std::vector<double> n_trans, cum_n_trans, h, w_sum;
    

    class invert_data{
      public:
        std::vector<int> in;
        std::vector<int> out;
    };

    std::vector<invert_data> lr_inv;

    double tol;
    int maxit;
    int it = 0;
    bool conv = false;
    bool inc_lik = false;
    double llike = R_NegInf;

    double calc_like();
    void calc_weight_sums();
    void run();


    // constructor
    length_bias(Rcpp::NumericVector lambda, Rcpp::IntegerVector l, 
                Rcpp::IntegerVector r, Rcpp::IntegerVector t, 
                Rcpp::IntegerVector R0, Rcpp::NumericVector del_t,
                double toler, int max_it) {

      maxit = max_it;
      tol = toler;
      n_int = lambda.length();
      n_obs = l.length();
      n = static_cast<double>(n_obs);

      left = Rcpp::as< std::vector<int> >(l);
      right = Rcpp::as< std::vector<int> >(r);
      trun = Rcpp::as<std::vector<int> >(t);
      h = Rcpp::as< std::vector<double> >(lambda);
      risk_0 = Rcpp::as< std::vector<double> >(R0);
      delta_t = Rcpp::as< std::vector<double> >(del_t);


      n_trans.resize(n_int);
      cum_n_trans.resize(n_int + 1);
      w_sum.resize(n_int);
      surv.resize(n_int + 1, 1);
      G.resize(n_int + 1);

      // initiate cum_lambda
      for (int j = 1; j < n_int + 1; j++){
          surv[j] = surv[j - 1] * (1 - h[j - 1]);
      }
      surv[n_int] = 0;
      h[n_int - 1] = 1;

      // survival integral
      for (int j = n_int - 1; j >=0; j--){
          G[j] = G[j + 1] + surv[j] * delta_t[j];
      }

      /* Invert Data */
      lr_inv.resize(n_int + 1);
      std::vector<int> left_sizes(n_int + 1);
      std::vector<int> right_sizes(n_int + 1);

      // find sizes
      for (int i = 0; i < n_obs; i++) {
        left_sizes[left[i]]++;
        right_sizes[right[i]]++;
      }

      // make data_inv correct size
      for (int j = 0; j < n_int + 1; j++) {
        lr_inv[j].in.resize(left_sizes[j]);
        lr_inv[j].out.resize(right_sizes[j]);
      }

      int curr_l, curr_r;
      int curr_j_l, curr_j_r;
      // transpose data
      for (int i = 0; i < n_obs; i++) {
        curr_l = left[i];
        curr_r = right[i];
        curr_j_l = left_sizes[curr_l] - 1;
        curr_j_r = right_sizes[curr_r] - 1;
        lr_inv[curr_l].in[curr_j_l] = i;
        lr_inv[curr_r].out[curr_j_r] = i;

        left_sizes[curr_l]--;
        right_sizes[curr_r]--;
      }
    };

};
