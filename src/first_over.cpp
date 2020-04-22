#include <Rcpp.h>
using namespace Rcpp;
#include "vector_sort.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' first_over function
//' Takes in a sorted NumericVector and a target double. Find the first element
//' in the vector which is larger than the target. Returns the index. Recall R 
//' starts indices at 1 and this is what the function uses. 
//' 
//' @param x: A sorted NumericVector
//' @param target: A double.
//' @return An integer for the index of the first element in x that is larger than the target.
//' @examples x<-runif(0,1,n=5) %>% sort() 
//' target=0.2 
//' find_index(x,target)
//' @export
// [[Rcpp::export]]
int first_over(NumericVector x, double target) {
  //check input
  if(is_sorted(x,true)==false){
    stop("The elements of the position vector are not in ascending order.");
  }
  
  NumericVector rel_dist = x - target;
  //Rcout<<rel_dist<<std::endl; 
  int n = x.size();
  for (int i = 0; i < n ; i++){
    if ( rel_dist[i] >= 0){
      return (i+1); //account for indices starting at 1 in R.
    }
  }
  warning("first_over. All elements of input vector are less than the target. Returning NA.");
  return R_NaN;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x = runif(5,0,1) %>% sort()
first_over(x,0.5)
*/
