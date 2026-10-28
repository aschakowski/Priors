data{
  
  // total number of nominated clusters
  int T;
  
  // number of alternative clusters (including the chosen one)
  int N_clusters;

} generated quantities {
   vector[T] log_lik;

    for(t in 1:T){
      
      // Likelihood evaluation //

      // Cluster nomination
      //Add log probability of observed block choice to target, note that chosen spot is last one in matrix
      log_lik[t] = categorical_logit_lpmf(N_clusters|rep_vector(0, N_clusters)); 
      
    }

}
