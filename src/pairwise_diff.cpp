#include <Rcpp.h>
#include "matrix_check.h"
using namespace Rcpp;

//' pairwise_diff 
//' 
//' Computes the number of pairwise differences between all the rows of a 
//' NumericMatrix. 
//' 
//' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
//' @return scalar value of the number of pairwise differences
//' @examples G = matrix(sample(0:1, size =5*5 , replace = TRUE), nc = 5), pairwise_diff(G) 
//' @export
// [[Rcpp::export]]
int pairwise_diff(NumericMatrix G) {
  
  //check input. 
  if(have_na(G)){
    Rcpp::warning("Input matrix is NA. Returning NA.");
    return (R_NaN);
  }
  
  int num_sam=G.nrow();
  int num_seg=G.ncol();
  
  int total_pairs=num_sam*(num_sam-1)/2;
  NumericVector diff(total_pairs);
  
  int k = 0;
  
  for (int n=0; n<=num_sam-1; n++){
    
    //loop through each pair
    for (int i=n+1; i<=num_sam-1;i++){
      
      //loop through each SNP
      for (int j=0; j<=num_seg-1; j++){
        
        if(G(n,j)!=G(i,j)){
          //index -1 because indices start at 0
          diff(k)+=1;
        }
      }
      //      Rcout<<"Pair"<<k<<" diff"<<diff(k)<<std::endl;
      k+=1;
    }
  }
  
  int ans = sum(diff);
  return ans;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
