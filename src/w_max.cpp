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

//' w_max function
//' 
//' Computes w_max for a binary, genome matrix. Each column is a SNP, each row is a sample. 
//' 
//' @param G: Binary genome matrix consisting of 1's and 0's.
//' @return scalar value of w_max for the genome matrix
//' @examples w_max(G)
//' @export
// [[Rcpp::export]]
NumericVector w_max(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
