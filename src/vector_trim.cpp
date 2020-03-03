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

//' vector_trim function
//' 
//' Trims outer elements of a vector based on some provided central index. 
//' 
//' @param x: NumericVector to trim
//' @param cen: Integer index for central element
//' @param k: Number of elements to retain from both sides of the center.
//' @return A NumericVector with outer elements removed. Up to k elements can be 
//' retained from both sides of the center. 
//' @examples x<-seq(1,5,by=1)
//' vector_trim(x,cen=3,k=1)
//' @export
// [[Rcpp::export]]
NumericVector vector_trim(NumericVector x, int cen, int k) {
  //c indices start at 0
  cen = cen-1;
  int last_index=x.length()-1; 
  
  int start = std::max(0 , cen-k);
  int end = std::min(last_index , cen+k); 
  NumericVector out = x[Rcpp::Range(start,end)];
  
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
