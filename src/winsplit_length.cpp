#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' winsplit_length function
//' 
//' Splits a genome matrix into equally sized windows as determined by the number of bases. 
//' The windows can differ in their number of columns when SNPs are not uniformly distributed
//' across the region. 
//' 
//' @param G: A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. Each column is a SNP. Each row is a sampled individual. 
//' @param pos: A NumericVector consisting of values between 0,1. The ith element is the position of the ith SNP in G. 
//' @param n: integer for the number of desired windows in output list. 
//' @return A NumericMatrix list of the windows
//' @examples  
//' seq <-matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
//' pos<-runif(0,1,n=5) %>% sort()
//' winsplit_length(seq,pos)
// [[Rcpp::export]]
List winsplit_length(NumericMatrix G,NumericVector pos, int n) {
  //compute genome length of each window
  double len=1/n;
  //find indices for the starting points of each window
  
  return 0;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
