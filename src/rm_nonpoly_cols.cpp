#include <Rcpp.h>
using namespace Rcpp;
#include "vec_equal.h"


//' nonpoly_cols function 
//' 
//' Takes a NumericMatrix and returns a list. The first element is a thinner matrix where
//' all the columns with the same element are removed. The second element is a vector of 
//' the indices of any columns that have all 0's, 1's and NA's. If none such columns exist, 
//' the vector is NA. 
//' 
//' @param G: A binary NumericMatrix of 0's and 1's. 
//' @return A list. Element 1 is a NumericMatrix with subset columns of G. Any columns with only 1's or only 0's
//' have been removed. Element 2 is vector containing the column indices of the polymorphic
//' columns. Note we use 1 indexing here for alignment with R.  
//' @examples G <-matrix(sample(0:1, size = 40, replace = TRUE), nc = 8)
//' rm_nonpoly_cols(G)
//' @export
// [[Rcpp::export]]

List rm_nonpoly_cols(NumericMatrix G) {
  
  int cols = G.ncol();
  int rows = G.nrow();
  
  //initialise an output matrix H with NAs.
  NumericMatrix H(rows,cols);
  std::fill(H.begin(),H.end(),R_NaN);
  
  //initialise vector to store the indices of all non-polymorphic sites
  NumericVector P(cols);
  std::fill(P.begin(), P.end(), R_NaN);
  
  //set a counter for the index of the next copied column
  int col_count = 0;

  for(int i=0; i < cols; i++){
    NumericVector Gcol = G(_,i);
    // Rcout << Gcol << std::endl; 
    // Rcout << NumericVector(rows,0.0) <<std::endl; 
    
    //check if the column has either all 0's or all 1's.
    bool all_0 = vec_equal(Gcol, NumericVector(rows,0.0));
    bool all_1 = vec_equal(Gcol, NumericVector(rows,1.0));
    bool all_na = all(is_na(Gcol));
    
    if(all_0 ==TRUE || all_1==TRUE || all_na == TRUE){
      continue;
    } else {
      H(_,col_count) = Gcol;
      P(col_count) = i;
      col_count = col_count + 1;
    }
  }
  
  //throw error if H has no columns
  if(col_count == 0){
    warning("Warning: All columns in genome matrix are non-polymorphic. Returning NA");
    List output = List::create(R_NaN,R_NaN);
    return output;
  }
  
  
  NumericMatrix mat = H( Range(0,rows - 1), Range(0, col_count - 1));
  // + 1 to accomodate 1 indexing in R
  NumericVector indices = P[Rcpp::Range(0, col_count - 1)] + 1; 
  
  List output = List::create(mat, indices);
  return output;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
