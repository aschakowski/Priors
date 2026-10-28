functions{
  
  // within chain parallelization
  real partial_sum_lpmf(array[] int slice,
                        int start,
                        int end,
                        array[] matrix Distances_social,
                        array[] matrix Distances_successes,
                        array[] matrix Distances_losses,
                        matrix Distances_ecological_knowledge,
                        matrix Locality,
                        vector Time_since_start,
                        array[] int success_columns,
                        array[] int loss_columns,
                        vector weights,
                        vector lambdas,
                        int N_spots,
                        array[] int lake,
                        array[] int id,
                        array[] int trip,
                        matrix v_lake,
                        matrix v_id,
                        matrix v_trip,
                        vector v_ones,
                        int columns) {
    
    real target_tmp=0;
    
    for(t in start:end){
      
      
      // Preallocation of feature vectors // 

      //preallocate social feature
      vector[N_spots] social_feature;
      
      // preallocate wins feature
      vector[N_spots] success_feature;
      
      // preallocate losses feature
      vector[N_spots] loss_feature;
      
      // preallocate locality feature
      vector[N_spots] locality_feature;
      
      // preallocate locality feature
      vector[N_spots] expectation_feature;
      
      // Computation of Feature values // 

      // compute social feature
      social_feature = (exp( - (Distances_social[t] .^ 2) / (2*square(exp(lambdas[1] + v_lake[lake[t], 9])))) .* (1.0 / (1 + exp(-20 * (Distances_social[t] - 5))))) * v_ones;
      
      // compute success feature
      if (Distances_successes[t][1, 1] < 0){
        success_feature = rep_vector(0, N_spots);
      } else {
        success_feature = exp( - (Distances_successes[t][, 1:success_columns[t]] .^ 2) /  (2*square(exp(lambdas[2] + v_lake[lake[t], 10])))) * rep_vector(1.0, success_columns[t]); //rep_vector(1.0/columns, loss_columns[t]);
      }
      
      // compute loss feature
      if (Distances_losses[t][1, 1] < 0){
        loss_feature = rep_vector(0, N_spots);
      } else {
        loss_feature = exp( - (Distances_losses[t][, 1:loss_columns[t]] .^ 2) /  (2*square(exp(lambdas[2] + v_lake[lake[t], 10])))) * rep_vector(1.0, loss_columns[t]); //rep_vector(1.0/columns, loss_columns[t]);
      }
      
      // compute locality feature
      locality_feature = Locality[, t] / 1000;
      
      // compute expectation feature
      expectation_feature = Distances_ecological_knowledge[, t] / 1000;
      
      // Likelihood evaluation //

      // Spot selection: different for first spot
      target_tmp += categorical_logit_lpmf(1|(weights[1] + v_trip[trip[t], 1] + v_lake[lake[t],1] + v_id[id[t], 1])*locality_feature +
                                             (weights[2] + v_trip[trip[t], 2] + v_lake[lake[t],2] + v_id[id[t], 2])*(locality_feature*Time_since_start[t]) +
                                             (weights[3] + v_trip[trip[t], 3] + v_lake[lake[t],3] + v_id[id[t], 3])*success_feature + 
                                             (weights[4] + v_trip[trip[t], 4] + v_lake[lake[t],4] + v_id[id[t], 4])*loss_feature + 
                                             (weights[5] + v_trip[trip[t], 5] + v_lake[lake[t],5] + v_id[id[t], 5])*social_feature + 
                                             (weights[6] + v_trip[trip[t], 6] + v_lake[lake[t],6] + v_id[id[t], 6])*expectation_feature + 
                                             (weights[7] + v_trip[trip[t], 7] + v_lake[lake[t],7] + v_id[id[t], 7])*expectation_feature .* success_feature + 
                                             (weights[8] + v_trip[trip[t], 8] + v_lake[lake[t],8] + v_id[id[t], 8])*expectation_feature .* loss_feature);
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
  
  // matrix containing distances to all other participants for each spot
  array[T] int success_columns; //how many columns to subset from distance matrix
  array[T] int loss_columns; //how many columns to subset from distance matrix
  
  int columns;
  array[T] matrix[N_spots, columns] Distances_social;
  array[T] matrix[N_spots, max(success_columns)] Distances_successes;
  array[T] matrix[N_spots, max(success_columns)] Distances_losses;
  
  // locality
  matrix[N_spots, T] Locality;
  
  //priors
  matrix[N_spots, T] Distances_ecological_knowledge;

  // time since start
  vector[T] Time_since_start;
  
  // choices for reduce sum
  array[T] int choices;
  
  //parallelization
  int grainsize;
  
} transformed data {
  
  // vector of ones for fast computation of rowMeans
  vector[columns] v_ones = rep_vector(1.0, columns); //rep_vector(1.0/columns, columns);
  
}

parameters{
  
  // Feature selection
  // bandwidth for feature selection
  vector[2] lambdas;
  
  // Spot selection
  // feature weights
  vector[8] weights;
  
  // varying effects: lakes
  matrix[10, N_lakes] z_lake;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[10] sigma_lake;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[10] Rho_lake; //Cholesky factor for covariance of learning parameters among individuals
  
  // varying effects: unique participants on lake
  matrix[8, N_ids] z_id;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[8] sigma_id;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[8] Rho_id; //Cholesky factor for covariance of learning parameters among individuals
  
  // varying effects: spots
  matrix[8, N_trips] z_trip;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[8] sigma_trip;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[8] Rho_trip; //Cholesky factor for covariance of learning parameters among individuals

} transformed parameters {
  
  // non centered parameterization
  matrix[N_lakes, 10] v_lake;
  v_lake = (diag_pre_multiply(sigma_lake, Rho_lake) * z_lake)';

  matrix[N_ids, 8] v_id;
  v_id = (diag_pre_multiply(sigma_id, Rho_id) * z_id)';
  
  matrix[N_trips, 8] v_trip;
  v_trip = (diag_pre_multiply(sigma_trip, Rho_trip) * z_trip)';
  
}

model{
  
  // Feature selection priors
  // priors on feature weights (normal(0,1) informative prior after exp)
  lambdas ~ normal(3, .5);
  
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
                       Distances_social,
                       Distances_successes,
                       Distances_losses,
                       Distances_ecological_knowledge,
                       Locality,
                       Time_since_start,
                       success_columns,
                       loss_columns,
                       weights,
                       lambdas,
                       N_spots,
                       lake,
                       id,
                       trip,
                       v_lake,
                       v_id,
                       v_trip,
                       v_ones,
                       columns);
  
} generated quantities {
   
  vector[T] log_lik;

    for(t in 1:T){
      
      // Preallocation of feature vectors // 

      //preallocate social feature
      vector[N_spots] social_feature;
      
      // preallocate wins feature
      vector[N_spots] success_feature;
      
      // preallocate losses feature
      vector[N_spots] loss_feature;
      
      // preallocate locality feature
      vector[N_spots] locality_feature;
      
      // preallocate locality feature
      vector[N_spots] expectation_feature;
      
      // Computation of Feature values // 

      // compute social feature
      social_feature = (exp( - (Distances_social[t] .^ 2) / (2*square(exp(lambdas[1] + v_lake[lake[t], 9])))) .* (1.0 / (1 + exp(-20 * (Distances_social[t] - 5))))) * v_ones;
      
      // compute success feature
      if (Distances_successes[t][1, 1] < 0){
        success_feature = rep_vector(0, N_spots);
      } else {
        success_feature = exp( - (Distances_successes[t][, 1:success_columns[t]] .^ 2) /  (2*square(exp(lambdas[2] + v_lake[lake[t], 10])))) * rep_vector(1.0, success_columns[t]); //rep_vector(1.0/columns, loss_columns[t]);
      }
      
      // compute loss feature
      if (Distances_losses[t][1, 1] < 0){
        loss_feature = rep_vector(0, N_spots);
      } else {
        loss_feature = exp( - (Distances_losses[t][, 1:loss_columns[t]] .^ 2) /  (2*square(exp(lambdas[2] + v_lake[lake[t], 10])))) * rep_vector(1.0, loss_columns[t]); //rep_vector(1.0/columns, loss_columns[t]);
      }
      
      // compute locality feature
      locality_feature = Locality[, t] / 1000;
      
      // compute expectation feature
      expectation_feature = Distances_ecological_knowledge[, t] / 1000;
      
      // Likelihood evaluation //

      // Spot selection: different for first spot
      log_lik[t] = categorical_logit_lpmf(1|(weights[1] + v_trip[trip[t], 1] + v_lake[lake[t],1] + v_id[id[t], 1])*locality_feature +
                                             (weights[2] + v_trip[trip[t], 2] + v_lake[lake[t],2] + v_id[id[t], 2])*(locality_feature*Time_since_start[t]) +
                                             (weights[3] + v_trip[trip[t], 3] + v_lake[lake[t],3] + v_id[id[t], 3])*success_feature + 
                                             (weights[4] + v_trip[trip[t], 4] + v_lake[lake[t],4] + v_id[id[t], 4])*loss_feature + 
                                             (weights[5] + v_trip[trip[t], 5] + v_lake[lake[t],5] + v_id[id[t], 5])*social_feature + 
                                             (weights[6] + v_trip[trip[t], 6] + v_lake[lake[t],6] + v_id[id[t], 6])*expectation_feature + 
                                             (weights[7] + v_trip[trip[t], 7] + v_lake[lake[t],7] + v_id[id[t], 7])*expectation_feature .* success_feature + 
                                             (weights[8] + v_trip[trip[t], 8] + v_lake[lake[t],8] + v_id[id[t], 8])*expectation_feature .* loss_feature);
    } 

}