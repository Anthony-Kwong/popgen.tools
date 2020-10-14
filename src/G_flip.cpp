#include <Rcpp.h>
#include "G_flip.h"
using namespace Rcpp;

//' G_flip function
//' 
//' Denoises a genome matrix by flipping an element in each column. For each column, 
//' we flip the element which results in the smallest number of pairwise differences
//' in the whole matrix. If there is a tie, we flip the element closest to the top. 
//' 
//' @param: G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
//' @return: A NumericMatrix. 
//' @export
// [[Rcpp::export]]
NumericMatrix G_flip(NumericMatrix G) {
  int snp = G.ncol();
  Rcpp::NumericMatrix K(Rcpp::clone(G));
  //NumericMatrix H = Gcol_denoise(K,2);
  for(int i=0;i<snp;i++){
    NumericMatrix H = Gcol_denoise(K,i);
    K = H;
  }
  return K;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
