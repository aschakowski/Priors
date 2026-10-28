functions{
  
  // within chain parallelization
  real partial_sum_lpmf(array[] int slice,
                        int start,
                        int end,
                        int N_spots) {
    
    real target_tmp=0;
    
    for(t in start:end){
      
      // Spot selection: random. remove first spots from fitting
      target_tmp += categorical_logit_lpmf(1|rep_vector(0, N_spots));
      
    }
  
    return target_tmp;

  }
  
} data{
  
  // length of complete time series (i.e., spot choices).
  int T;
  
  // number of available angling spots (including the chosen one)
  int N_spots;
  
  // choices for reduce sum
  array[T] int choices;

  //parallelization
  int grainsize;
  
} 

parameters{
  
}

model{
  
  
  // add to ll
  target += reduce_sum(partial_sum_lpmf,
                       choices,
                       grainsize,
                       N_spots);
  
} generated quantities {
   
  vector[T] log_lik;
  
    for(t in 1:T){
      
      log_lik[t] = categorical_logit_lpmf(1|rep_vector(0, N_spots));
      
    }
  
}