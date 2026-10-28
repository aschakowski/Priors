functions{
  
  // within chain parallelization
  real partial_sum_lpmf(array[] int slice,
                        int start,
                        int end,
                        matrix Locality,
                        vector Time_since_start,
                        vector weights,
                        int N_spots,
                        array[] int lake,
                        array[] int id,
                        array[] int trip,
                        matrix v_lake,
                        matrix v_id,
                        matrix v_trip) {
    
    real target_tmp=0;
    
    for(t in start:end){
      
      
      // Preallocation of feature vectors // 

      // preallocate locality feature
      vector[N_spots] locality_feature;
      
      // Computation of Feature values // 

      // compute locality feature
      locality_feature = Locality[, t] / 1000;
      
      // Likelihood evaluation //

      // Spot selection: different for first spot
      target_tmp += categorical_logit_lpmf(1|(weights[1] + v_trip[trip[t], 1] + v_lake[lake[t],1] + v_id[id[t], 1])*locality_feature + 
                                             (weights[2] + v_trip[trip[t], 2] + v_lake[lake[t],2] + v_id[id[t], 2])*(locality_feature*Time_since_start[t]));

    }
  
    return target_tmp;

  }
  
} data{
  
  // length of complete time series (i.e., spot choices).
  int T;
  
  // number of available angling spots (including the chosen one)
  int N_spots;
  
  // number of foraging trips and index
  array[T] int trip;
  int N_trips;
  
  // number of participants and index
  array[T] int id;
  int N_ids;
  
  // number of lakes and index
  array[T] int lake;
  int N_lakes;
  
  // locality
  matrix[N_spots, T] Locality;
  
  // time since start
  vector[T] Time_since_start;
  
  // choices for reduce sum
  array[T] int choices;
  
  //parallelization
  int grainsize;
    
} 

parameters{
  
  // Spot selection
  // feature weights
  vector[2] weights;
  
  // varying effects: lakes
  matrix[2, N_lakes] z_lake;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[2] sigma_lake;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[2] Rho_lake; //Cholesky factor for covariance of learning parameters among individuals
  
  // varying effects: unique participants on lake
  matrix[2, N_ids] z_id;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[2] sigma_id;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[2] Rho_id; //Cholesky factor for covariance of learning parameters among individuals
  
  // varying effects: spots
  matrix[2, N_trips] z_trip;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[2] sigma_trip;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[2] Rho_trip; //Cholesky factor for covariance of learning parameters among individuals

} transformed parameters {
  
  // non centered parameterization
  matrix[N_lakes, 2] v_lake;
  v_lake = (diag_pre_multiply(sigma_lake, Rho_lake) * z_lake)';

  matrix[N_ids, 2] v_id;
  v_id = (diag_pre_multiply(sigma_id, Rho_id) * z_id)';
  
  matrix[N_trips, 2] v_trip;
  v_trip = (diag_pre_multiply(sigma_trip, Rho_trip) * z_trip)';
  
}

model{
  
  // spot selection priors
  //Define prior distribution of varying individual effects
  to_vector(z_lake) ~ normal(0, 1);
  sigma_lake ~ exponential(1);
  Rho_lake ~ lkj_corr_cholesky(4);
    
  //Define prior distribution of varying individual effects
  to_vector(z_id) ~ normal(0, 1);
  sigma_id ~ exponential(1);
  Rho_id ~ lkj_corr_cholesky(4);

  //Define prior distribution of varying individual effects
  to_vector(z_trip) ~ normal(0, 1);
  sigma_trip ~ exponential(1);
  Rho_trip ~ lkj_corr_cholesky(4);
    
  // priors on feature weights (normal(0,1) informative prior after softmax)
  //to_vector(weights) ~ normal(0,1);
  to_vector(weights) ~ normal(0,5);
  
  // add to ll
  target += reduce_sum(partial_sum_lpmf,
                       choices,
                       grainsize,
                       Locality,
                       Time_since_start,
                       weights,
                       N_spots,
                       lake,
                       id,
                       trip,
                       v_lake,
                       v_id,
                       v_trip);
  
} generated quantities {
   
  vector[T] log_lik;

    for(t in 1:T){
      
      
      // preallocate locality feature
      vector[N_spots] locality_feature;
      
      // compute locality feature
      locality_feature = Locality[, t] / 1000;
      
      // Likelihood evaluation //

      // Spot selection: different for first spot
      log_lik[t] = categorical_logit_lpmf(1|(weights[1] + v_trip[trip[t], 1] + v_lake[lake[t],1] + v_id[id[t], 1])*locality_feature + 
                                            (weights[2] + v_trip[trip[t], 2] + v_lake[lake[t],2] + v_id[id[t], 2])*(locality_feature*Time_since_start[t]));
    }
  
}