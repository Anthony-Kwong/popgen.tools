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

//' is_sorted function
//' 
//' Checks if a NumericVector is sorted numerically. 
//' 
//' @param x: NumericVector
//' @param ascend: A boolean. If TRUE, checks if the elements of x are in ascending order.
//' If FALSE, checks if they are in descending order. 
//' @return A boolean indicating whether the vector is sorted.
//' @examples x = c(1,3,5,7)
//' is_sorted(x, ascend = T)
//' @export 
// [[Rcpp::export]]
bool is_sorted(NumericVector x, bool ascend) {
  int n= x.size();
  if (ascend==true){
    for (int i = 0; i < (n-1); i++  ){
      if( x[i] > x[i+1]) {
        return false;
      }
    }
    return true;
  } else {
    for(int i = 0; i < (n-1); i++){
      if(x[i] < x[i+1]) {
        return false;
      }
    }
    return true;
  }
  stop("is_sorted: Check loop busted. This line should never be called. :( ");
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
