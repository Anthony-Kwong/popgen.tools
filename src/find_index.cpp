#include <Rcpp.h>
using namespace Rcpp;

//' find_index function
//' 
//' Searches a NumericVector for the element that is closest to a target value. 
//' Returns the corresponding index. 
//' @param x: A NumericVector with the elements are arranged in ascending order. 
//' @param target: a double
//' @return An integer for the index of the closest element in NumericVector x
//' @examples
//' x<-runif(0,1,n=5) %>% sort()
//' target=0.2

// [[Rcpp::export]]
int find_index(NumericVector x,double target) {
  NumericVector y=x-target;
  //distances must be non-negative
  NumericVector dist=Rcpp::abs(y);
  int index=which_min(dist);
  return index;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(1350)
x<-runif(0,1,n=5) %>% sort()
d<-0.2
find_index(x,d)
*/
