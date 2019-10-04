#This function was written to test the theta_t function. 

#Count the number of pair wise differences between rows of a numeric matrix, normalised 
#by the number of pairs available. 

#input: Binary matrix G consisting of 1's and 0's. 
#output: Number of pairwise differences normalised by the number of pairs. 

test_theta_t<-function(G){
  
  pair_diff=dist(G,method="manhattan")
  total_diff=sum(pair_diff)
  nsam=nrow(G)
  total_pairs=choose(nsam,2)
  
  return(total_diff/total_pairs)
}
