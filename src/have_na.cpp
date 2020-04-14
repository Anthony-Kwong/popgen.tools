#include <Rcpp.h>
using namespace Rcpp;

//' have_na function
//' 
//' Checks if a NumericMatrix has any NA values. Returns T if it does, F otherwise.
//' 
//' @param M: A NumericMatrix
//' @return A boolean indicating whether the input matrix has any NAs. 
//' @examples A = matrix(1:4, nrow = 2)
//' have_na(A)
//' @export 
// [[Rcpp::export]]
bool have_na (NumericMatrix M) {
  return is_true(any(is_na(M)));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
A = matrix(1:4, nrow = 2)
have_na(A)
*/
