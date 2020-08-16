#include <Rcpp.h>
using namespace Rcpp;

//' matrix_subset function
//' 
//' Subsets selected columns from a NumericMatrix. 
//' 
//' @param G: A NumericMatrix
//' @param y: An IntegerVector of indices. Uses 1 indexing to suit R.
//' @return A NumericMatrix with subset columns of G using the indices in y.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_subset(Rcpp::NumericMatrix G, Rcpp::IntegerVector y) { 
  
  // Determine the number of observations
  int n_cols_out = y.size();
  
  // Create an output matrix
  Rcpp::NumericMatrix out = Rcpp::no_init(G.nrow(), n_cols_out);
  
  // Loop through each column and copy the data. 
  for(unsigned int z = 0; z < n_cols_out; z++) {
    out(Rcpp::_, z) = G(Rcpp::_, y[z]-1); //-1 to change to 1 indexing in R
  }
  return out;
}

//adapted from https://stackoverflow.com/questions/62118084/rcpp-select-subset-numericmatrix-column-by-a-numericvector

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
