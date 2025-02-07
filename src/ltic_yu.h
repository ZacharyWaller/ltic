class ltic_yu{
  public:
    int n_int;
    int n_obs;
    double n;


    std::vector<int> left;
    std::vector<int> right;
    std::vector<int> trun;
    std::vector<double> s_0;
    std::vector<double> s_1;

    std::vector<double> deriv_1;
    std::vector<double> deriv_2;
    std::vector<double> cum_lambda;
    std::vector<double> surv;
    std::vector<double> n_trans, cum_n_trans, h, w_sum;

    class invert_data{
      public:
        std::vector<int> in;
        std::vector<int> out;
        std::vector<int> trun_time;
    };

    std::vector<invert_data> lr_inv;

    double tol;
    int maxit;
    int it = 0;
    bool conv = false;
    double llike = R_NegInf;

    double calc_like();
    void em_algo();
    void calc_weight_sums();
    void convert_to_haz();
    void convert_to_surv();
    void newton_algo();
    void calc_derivs();
    void half_steps();
    void run();


    // constructor
    ltic_yu(Rcpp::NumericVector lambda, Rcpp::IntegerVector l, 
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
      h = Rcpp::as<std::vector<double> > (lambda);
      s_0 = Rcpp::as< std::vector<double> >(lambda);
      s_1 = s_0;

      deriv_1.resize(n_int + 1);
      deriv_2.resize(n_int + 1);
      cum_lambda.resize(n_int + 1);
      surv.resize(n_int + 1, 1);
      n_trans.resize(n_int);
      cum_n_trans.resize(n_int + 1);
      w_sum.resize(n_int);

      /* Initiate Cumulative Hazards */
      // for (int j = 1; j < n_int + 1; j++){
      //     cum_lambda[j] = cum_lambda[j - 1] + lambda[j - 1];
      // }
      cum_lambda[0] = 0;
      cum_lambda[n_int] = R_PosInf;

      // initiate surv
      surv[0] = 1;
      for (int j = 0; j < n_int; j++){
          surv[j + 1] = surv[j] - s_0[j];
      }

      /* Invert Data */
      lr_inv.resize(n_int + 1);
      std::vector<int> left_sizes(n_int + 1);
      std::vector<int> right_sizes(n_int + 1);
      std::vector<int> trun_sizes(n_int + 1);

      // find sizes
      for (int i = 0; i < n_obs; i++) {
        left_sizes[left[i]]++;
        right_sizes[right[i]]++;
        trun_sizes[trun[i]]++;
      }

      // make data_inv correct size
      for (int j = 0; j < n_int + 1; j++) {
        lr_inv[j].in.resize(left_sizes[j]);
        lr_inv[j].out.resize(right_sizes[j]);
        lr_inv[j].trun_time.resize(trun_sizes[j]);
      }

      int curr_l, curr_r, curr_t;
      int curr_j_l, curr_j_r, curr_j_t;
      // transpose data
      for (int i = 0; i < n_obs; i++) {
        curr_l = left[i];
        curr_r = right[i];
        curr_t = trun[i];
        curr_j_l = left_sizes[curr_l] - 1;
        curr_j_r = right_sizes[curr_r] - 1;
        curr_j_t = trun_sizes[curr_t] - 1;
        lr_inv[curr_l].in[curr_j_l] = i;
        lr_inv[curr_r].out[curr_j_r] = i;
        lr_inv[curr_t].trun_time[curr_j_t] = i;

        left_sizes[curr_l]--;
        right_sizes[curr_r]--;
        trun_sizes[curr_t]--;
      }

    };

};
