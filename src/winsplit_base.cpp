#include <Rcpp.h>
using namespace Rcpp;
#include "winsplit_length.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' winsplit_base function
//' 
//' Split a genome matrix into equally sized windows as determined by
//' the number of bases. The windows can differ in the number of columns
//' because SNPs are not always uniform across a given region. 
//' 
//' @param G: A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. Each column is a SNP. Each row is a sampled individual. 
//' @param pos: A NumericVector consisting of values between 0,1. The ith element is the position of the ith SNP in G.
//' @param n: integer for the number of desired windows in the output list. 
//' @return A NumericMatrix list of the windows
//' @examples seq <-matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
//' pos<-runif(0,1,n=5) %>% sort()
//' winsplit_length(seq,pos,5)
//' @export
// [[Rcpp::export]]
NumericVector winsplit_base(NumericMatrix G, NumericVector pos, int n) {
  //check inputs
  if(any_sug(pos<0) || any_sug(pos>1)){
    stop("Elements of pos must be between 0,1.");
  }
  
  Rcout<<find_index(pos,0.2)<<std::endl;
  return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
seq <-matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
pos<-c(0,0.19,0.200000001,0.5,0.9,2)
winsplit_base(seq,pos,n=3)
*/
