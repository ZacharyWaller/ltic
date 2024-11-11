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

    std::vector<double> surv;
    std::vector<double> w_sum;

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
    bool inc_lik = false;
    double llike = R_NegInf;

    double calc_like();
    void calc_weight_sums();
    void run();


    // constructor
    yu(Rcpp::NumericVector s, Rcpp::IntegerVector l, Rcpp::IntegerVector r, 
       Rcpp::IntegerVector t, double toler, int max_it) {

      maxit = max_it;
      tol = toler;
      n_int = s.length();
      n_obs = l.length();
      n = static_cast<double>(n_obs);

      left = Rcpp::as< std::vector<int> >(l);
      right = Rcpp::as< std::vector<int> >(r);
      trun = Rcpp::as<std::vector<int> >(t);
      s_0 = Rcpp::as< std::vector<double> >(s);
      s_1 = s_0;

      w_sum.resize(n_int);
      surv.resize(n_int + 1);

      // initiate surv
      surv[0] = 1;
      for (int j = 0; j < n_int; j++){
          surv[j + 1] = surv[j] - s[j];
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
