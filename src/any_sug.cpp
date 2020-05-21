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

//' any_sug function
//' 
//' A wrapper for the sugar function any. It takes a LogicalVector and returns T if any of the elements are T.
//' Source: https://gallery.rcpp.org/articles/sugar-any/
//' 
//' @param x: LogicalVector.
//' @return A boolean.
//' @examples  x <- c(3, 9, 0, 2, 7, -1, 6)
//' any_sug(x>9)
//' @export
// [[Rcpp::export]]
bool any_sug(LogicalVector x){
  // Note the use of is_true to return a bool type 
  return is_true(any(x == TRUE));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x<-c(0,1,3,5,7)
any_sug(x<2)
any_sug(x>8)
*/
