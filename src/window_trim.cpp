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

//' window_trim function
//' 
//' Trims the outer columns of a matrix. 
//' 
//' @param G: NumericMatrix
//' @param cen: An integer designating a column index. This will be the center column of the ouput matrix.
//' @param k: An integer for the number of columns to include from either side of the center column. 
//' @return A NumericMatrix with the outer columns of G removed. 
//' @examples window_trim(G,1,5)
//' @export
// [[Rcpp::export]]
NumericMatrix window_trim(NumericMatrix G,int cen,int k) {
  
  //add error if cen is not a valid index
  
  //account for indices starting at 0 in c++
  cen=cen-1;
  
  //take indices of the k columns of either side of the center
  int start=cen-k;
  int end=cen+k;
  
  if(start<0){
    stop("window_trim: Trying to include too many columns to the left of center.");
  }
  
  int cols=G.ncol();
  
  if(end>(cols-1)){
    stop("window_trim: Trying to include too many columns to the right of center.");
  }
  
  NumericMatrix G2=G(_,Range(start,end));
  return G2;
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
