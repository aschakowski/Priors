functions{
  
  // within chain parallelization
  real partial_sum_lpmf(array[] int slice,
                        int start, //trip id start
                        int end, //trip id end
                        array[] int trip_start, // row index where each trip starts
                        array[] int trip_end, // row index where each trip ends
                        array[] int state, //either 0 (stay) or 1 (leave)
                        array[] int catch, // either 0 (no catch) or 1 (catch)
                        vector cumulative_catch, //cumulative number of caught fish at each spot
                        vector success_feature, // spatial success feature (posterior prediction of spot selection model)
                        vector loss_feature, // spatial loss feature (posterior prediction of spot selection model)
                        vector knowledge_feature, // spatial roughness feature
                        vector social_feature, // spatial social feature
                        vector time_since_event, // time since start or fish catch (if cumulative_catch > 0)
                        vector angling_time, // total angling time at spot
                        array[] int confidence, // confidence estimate of spot
                        vector beta, // regression weights
                        matrix v_trip, // trip offsets
                        matrix v_lake, // lake offsets
                        matrix v_id, // id offsets
                        array[] int lakes, // lake identifier
                        array[] int ids // id identifier
                        ) {
    
    real target_tmp=0;
    
    // loop through all trips
    for (unique_trip in start:end){
        
        // loop through each timestamp of angling trip
        for (t in (trip_start[unique_trip]):trip_end[unique_trip]){
        
          if (confidence[t]>-99){
          // add ll
              if (cumulative_catch[t] == 0){
                  target_tmp += bernoulli_logit_lpmf( state[t] |  beta[1]+v_trip[unique_trip, 1]+v_lake[lakes[t], 1] + v_id[ids[t], 1] + 
                                                                  (beta[3]+v_trip[unique_trip, 3]+v_lake[lakes[t], 3] + v_id[ids[t], 3])*time_since_event[t] + 
                                                                  (beta[4]+v_trip[unique_trip, 4]+v_lake[lakes[t], 4] + v_id[ids[t], 4])*success_feature[t] + 
                                                                  (beta[5]+v_trip[unique_trip, 5]+v_lake[lakes[t], 5] + v_id[ids[t], 5])*loss_feature[t] + 
                                                                  (beta[6]+v_trip[unique_trip, 6]+v_lake[lakes[t], 6] + v_id[ids[t], 6])*social_feature[t] + 
                                                                  (beta[6 + confidence[t]]+v_trip[unique_trip, 6 + confidence[t]]+v_lake[lakes[t], 6 + confidence[t]] + v_id[ids[t], 6 + confidence[t]])*knowledge_feature[t]);
              } else {
                  target_tmp += bernoulli_logit_lpmf( state[t] |  beta[1]+v_trip[unique_trip, 1]+v_lake[lakes[t], 1] + v_id[ids[t], 1] +
                                                                  beta[2]+v_trip[unique_trip, 2]+v_lake[lakes[t], 2] + v_id[ids[t], 2] + 
                                                                  (beta[3]+v_trip[unique_trip, 3]+v_lake[lakes[t], 3] + v_id[ids[t], 3])*time_since_event[t] + 
                                                                  (beta[4]+v_trip[unique_trip, 4]+v_lake[lakes[t], 4] + v_id[ids[t], 4])*success_feature[t] + 
                                                                  (beta[5]+v_trip[unique_trip, 5]+v_lake[lakes[t], 5] + v_id[ids[t], 5])*loss_feature[t] + 
                                                                  (beta[6]+v_trip[unique_trip, 6]+v_lake[lakes[t], 6] + v_id[ids[t], 6])*social_feature[t] + 
                                                                  (beta[6 + confidence[t]]+v_trip[unique_trip, 6 + confidence[t]]+v_lake[lakes[t], 6 + confidence[t]] + v_id[ids[t], 6 + confidence[t]])*knowledge_feature[t]);
              }

              // count up
              //counter += 1;
          }
        }

    } // trip
  
    return target_tmp;
      
  } // function

}

data{
  
  int<lower = 1> T; // length of time series
  array[T] int<lower = 0, upper = 1>  state;
  array[T] int<lower = 0, upper = 1>  catch;
  
  // time varying predictors
  vector[T] time_since_event; 
  vector[T] cumulative_catch; 
  vector[T] angling_time; 
  
  // spatial features
  vector[T] success_feature; 
  vector[T] loss_feature; 
  vector[T] knowledge_feature; 
  vector[T] social_feature; 
  array[T] int confidence;
  
  // grainsize for within chain parallelization
  int grainsize;
  
  // number of foraging trips
  int N_trips;
  int N_lakes;
  int N_ids;
  
  // trip indices for wihtin chain parallelization
  array[N_trips] int trips;
  array[T] int lakes;
  array[T] int ids;
   // row indices to subset trips
  array[N_trips] int trip_start; // contains first index of trip
  array[N_trips] int trip_end; // contains last index of trip

}
parameters {

  // regression parameter
  vector[9] beta;

  // trip level offsets
  matrix[9, N_trips] z_trip;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[9] sigma_trip;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[9] Rho_trip; //Cholesky factor for covariance of learning parameters among individuals
  
  // lake level offsets
  matrix[9, N_lakes] z_lake;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[9] sigma_lake;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[9] Rho_lake; //Cholesky factor for covariance of learning parameters among individuals
  
  // id level offsets
  matrix[9, N_ids] z_id;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[9] sigma_id;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[9] Rho_id; //Cholesky factor for covariance of learning parameters among individuals
  
  
} transformed parameters {

  // non centered parameterization
  matrix[N_trips, 9] v_trip;
  v_trip = (diag_pre_multiply(sigma_trip, Rho_trip) * z_trip)';

  // non centered parameterization
  matrix[N_lakes, 9] v_lake;
  v_lake = (diag_pre_multiply(sigma_lake, Rho_lake) * z_lake)';

  // non centered parameterization
  matrix[N_ids, 9] v_id;
  v_id = (diag_pre_multiply(sigma_id, Rho_id) * z_id)';

} model {

  // trip level offsets
  to_vector(z_trip) ~ normal(0, 1);
  sigma_trip ~ exponential(1);
  Rho_trip ~ lkj_corr_cholesky(4);

  // trip level offsets
  to_vector(z_lake) ~ normal(0, 1);
  sigma_lake ~ exponential(1);
  Rho_lake ~ lkj_corr_cholesky(4);

  // trip level offsets
  to_vector(z_id) ~ normal(0, 1);
  sigma_id ~ exponential(1);
  Rho_id ~ lkj_corr_cholesky(4);

  // regression priors
  beta ~ normal(0, 3);
  
  // updating priors
  //alpha ~ normal(0, 1);

  // add to ll
  target += reduce_sum(partial_sum_lpmf,
                       trips,
                       grainsize,
                       trip_start,
                       trip_end,
                       state,
                       catch,
                       cumulative_catch,
                       success_feature,
                       loss_feature,
                       knowledge_feature,
                       social_feature,
                       time_since_event,
                       angling_time,
                       confidence,
                       beta,
                       v_trip,
                       v_lake,
                       v_id,
                       lakes,
                       ids);
                        
} generated quantities {
   
   vector[T] log_lik;

    // loop through all trips
    for (unique_trip in 1:N_trips){
        
        // loop through each timestamp of angling trip
        for (t in (trip_start[unique_trip]):trip_end[unique_trip]){
        
          if (confidence[t]>-99){
          // add ll
              if (cumulative_catch[t] == 0){
                  log_lik[t] = bernoulli_logit_lpmf( state[t] |  beta[1]+v_trip[unique_trip, 1]+v_lake[lakes[t], 1] + v_id[ids[t], 1] + 
                                                                  (beta[3]+v_trip[unique_trip, 3]+v_lake[lakes[t], 3] + v_id[ids[t], 3])*time_since_event[t] + 
                                                                  (beta[4]+v_trip[unique_trip, 4]+v_lake[lakes[t], 4] + v_id[ids[t], 4])*success_feature[t] + 
                                                                  (beta[5]+v_trip[unique_trip, 5]+v_lake[lakes[t], 5] + v_id[ids[t], 5])*loss_feature[t] + 
                                                                  (beta[6]+v_trip[unique_trip, 6]+v_lake[lakes[t], 6] + v_id[ids[t], 6])*social_feature[t] + 
                                                                  (beta[6 + confidence[t]]+v_trip[unique_trip, 6 + confidence[t]]+v_lake[lakes[t], 6 + confidence[t]] + v_id[ids[t], 6 + confidence[t]])*knowledge_feature[t]);
              } else {
                  log_lik[t] = bernoulli_logit_lpmf( state[t] |  beta[1]+v_trip[unique_trip, 1]+v_lake[lakes[t], 1] + v_id[ids[t], 1] +
                                                                  beta[2]+v_trip[unique_trip, 2]+v_lake[lakes[t], 2] + v_id[ids[t], 2] + 
                                                                  (beta[3]+v_trip[unique_trip, 3]+v_lake[lakes[t], 3] + v_id[ids[t], 3])*time_since_event[t] + 
                                                                  (beta[4]+v_trip[unique_trip, 4]+v_lake[lakes[t], 4] + v_id[ids[t], 4])*success_feature[t] + 
                                                                  (beta[5]+v_trip[unique_trip, 5]+v_lake[lakes[t], 5] + v_id[ids[t], 5])*loss_feature[t] + 
                                                                  (beta[6]+v_trip[unique_trip, 6]+v_lake[lakes[t], 6] + v_id[ids[t], 6])*social_feature[t] + 
                                                                  (beta[6 + confidence[t]]+v_trip[unique_trip, 6 + confidence[t]]+v_lake[lakes[t], 6 + confidence[t]] + v_id[ids[t], 6 + confidence[t]])*knowledge_feature[t]);
              }

              // count up
              //counter += 1;
          }
        }

    } // trip

}
