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
                        vector time_since_event, // time since start or fish catch (if cumulative_catch > 0)
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
            
            // add ll
            target_tmp += bernoulli_logit_lpmf( state[t] |  beta[1]+v_trip[unique_trip, 1]+v_lake[lakes[t], 1] + v_id[ids[t], 1]);

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
  vector[1] beta;

  // trip level offsets
  matrix[1, N_trips] z_trip;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[1] sigma_trip;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[1] Rho_trip; //Cholesky factor for covariance of learning parameters among individuals
  
  // lake level offsets
  matrix[1, N_lakes] z_lake;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[1] sigma_lake;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[1] Rho_lake; //Cholesky factor for covariance of learning parameters among individuals
  
  // id level offsets
  matrix[1, N_ids] z_id;           //Matrix for our latent individual samples (z scores)
  vector<lower = 0>[1] sigma_id;  //Standard deviation of learning parameters among individuals
  cholesky_factor_corr[1] Rho_id; //Cholesky factor for covariance of learning parameters among individuals
  
  
} transformed parameters {

  // non centered parameterization
  matrix[N_trips, 1] v_trip;
  v_trip = (diag_pre_multiply(sigma_trip, Rho_trip) * z_trip)';

  // non centered parameterization
  matrix[N_lakes, 1] v_lake;
  v_lake = (diag_pre_multiply(sigma_lake, Rho_lake) * z_lake)';

  // non centered parameterization
  matrix[N_ids, 1] v_id;
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
  
  // add to ll
  target += reduce_sum(partial_sum_lpmf,
                       trips,
                       grainsize,
                       trip_start,
                       trip_end,
                       state,
                       catch,
                       cumulative_catch,
                       time_since_event,
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
            
            // add ll
            log_lik[t] = bernoulli_logit_lpmf( state[t] |  beta[1]+v_trip[unique_trip, 1]+v_lake[lakes[t], 1] + v_id[ids[t], 1]);

        }

    }

}
