functions{
  
  // within chain parallelization
  real partial_sum_lpmf(array[] int slice,
                        int start,
                        int end,
                        int N_clusters,
                        vector weights,
                        matrix roughness,
                        matrix displacement,
                        matrix depth,
                        matrix boundary,
                        matrix shore,
                        array[] int lake,
                        array[] int trip,
                        array[] int id,
                        matrix v_trip,
                        matrix v_lake,
                        matrix v_id,
                        array[] int confidence) {
    
    real target_tmp=0;
    
    for(t in start:end){
      
      // Likelihood evaluation //

      // Cluster nomination
      //Add log probability of observed block choice to target, note that chosen spot is last one in matrix
      target_tmp += categorical_logit_lpmf(N_clusters|(weights[1] + v_trip[trip[t],1] + v_lake[lake[t],1] + v_id[id[t], 1])*displacement[,t]); // + v_trip[trip[t],3] + v_lake[lake[t],3] + v_id[id[t], 3]
      
    }
  
    return target_tmp;

  }
  
} data{
  
  // total number of nominated clusters
  int T;
  
  // number of alternative clusters (including the chosen one)
  int N_clusters;
  
  // number of foraging trips and index
  array[T] int trip;
  int N_trips;
  
  // number of participants and index
  array[T] int id;
  int N_ids;
  
  // number of lakes and index
  array[T] int lake;
  int N_lakes;
  
  // features
  matrix[N_clusters, T] roughness;
  matrix[N_clusters, T] displacement;
  matrix[N_clusters, T] depth;
  matrix[N_clusters, T] boundary;
  matrix[N_clusters, T] shore;
  
  // choices for reduce sum
  array[T] int choices;
  
  
  array[T] int confidence;
  
  //parallelization
  int grainsize;
  
} parameters{
  
  // Cluster nomination
  // feature weights
  vector[1] weights;
  
  // varying effects: foraging trip
  matrix[1, N_trips] z_trip;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[1] sigma_trip;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[1] Rho_trip; //Cholesky factor for covariance of learning parameters among individuals
  
  // varying effects: lakes
  matrix[1, N_lakes] z_lake;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[1] sigma_lake;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[1] Rho_lake; //Cholesky factor for covariance of learning parameters among individuals
  
  // varying effects: unique participants on lake
  matrix[1, N_ids] z_id;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[1] sigma_id;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[1] Rho_id; //Cholesky factor for covariance of learning parameters among individuals
  
} transformed parameters {
  
  // non centered parameterization
  matrix[N_trips, 1] v_trip;
  v_trip = (diag_pre_multiply(sigma_trip, Rho_trip) * z_trip)';

  matrix[N_lakes, 1] v_lake;
  v_lake = (diag_pre_multiply(sigma_lake, Rho_lake) * z_lake)';

  matrix[N_ids, 1] v_id;
  v_id = (diag_pre_multiply(sigma_id, Rho_id) * z_id)';

}

model{
  
  // cluster nomination priors
  //Define prior distribution of varying individual effects
  to_vector(z_trip) ~ normal(0, 1);
  sigma_trip ~ exponential(1);
  Rho_trip ~ lkj_corr_cholesky(4);

  //Define prior distribution of varying individual effects
  to_vector(z_lake) ~ normal(0, 1);
  sigma_lake ~ exponential(1);
  Rho_lake ~ lkj_corr_cholesky(4);
    
  //Define prior distribution of varying individual effects
  to_vector(z_id) ~ normal(0, 1);
  sigma_id ~ exponential(1);
  Rho_id ~ lkj_corr_cholesky(4);

  // priors on feature weights (normal(0,1) informative prior after softmax)
  to_vector(weights) ~ normal(0,3);
  
  // add to ll
  target += reduce_sum(partial_sum_lpmf,
                       choices,
                       grainsize,
                       N_clusters,
                       weights,
                       roughness,
                       displacement,
                       depth,
                       boundary,
                       shore,
                       lake,
                       trip,
                       id,
                       v_trip,
                       v_lake,
                       v_id,
                       confidence);
  
} generated quantities {
   vector[T] log_lik;

    for(t in 1:max(choices)){
      
      // Likelihood evaluation //

      // Cluster nomination
      //Add log probability of observed block choice to target, note that chosen spot is last one in matrix
      log_lik[t] = categorical_logit_lpmf(N_clusters|(weights[1] + v_trip[trip[t],1] + v_lake[lake[t],1] + v_id[id[t], 1])*displacement[,t]); // + v_trip[trip[t],3] + v_lake[lake[t],3] + v_id[id[t], 3]
      
    }

}
