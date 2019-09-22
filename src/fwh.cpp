#include <Rcpp.h>
using namespace Rcpp;
#include "popgen.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' theta_h function
//' Computes theta_h for a genome matrix
//' @param G: Binary genome matrix with 0's and 1's. Each column is a SNP, each row is an individual.
//' @return scalar value of theta_h for that sampled population
//' @examples
//' @export
// [[Rcpp::export]]
double theta_h(NumericMatrix G){
  return 0;
}

//' fwh function
//' Computes Fay and Wu's H for a genome matrix. 
//' @param G: Binary genome matrix with 0's and 1's. Each column is a SNP, each row is an individual.
//' @return scalar value of Fay and Wu's H for that sampled population.  
//' @examples
//' @export
// [[Rcpp::export]]
NumericVector fwh(NumericVector x) {
  return x * 2;
}

//H is the difference between theta_pi (average het) and theta_h

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
