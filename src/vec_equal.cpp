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
//' vec_equal function for cpp
//' @param x,y: two numeric vectors of the same dimensions
//' @return boolean, true if the two vectors are the same. False otherwise. 
//' @examples vec_equal(x,y)
//' @export
// [[Rcpp::export]]
bool vec_equal(NumericVector x, NumericVector y){
  int size_x=x.size();
  int size_y=y.size();
  
  if(size_x != size_y){
    try{
      throw 1;
    }
    catch (int e){
      Rcout<<"Exception error # "<<e<<std::endl;
      Rcout<<"Vectors must have the same dimension!"<<std::endl;
      return FALSE;
    }
  }

  
  for(int i=0; i<size_x; i++){
    if(x[i]!=y[i]){
      return FALSE;
    }
  }
//  Rcout<<size_x<<" "<<size_y<<std::endl;
  return (TRUE);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x<-c(1,3,4)
y<-c(1,2,3)
z<-c(1,2,3,4)
vec_equal(y,z)
*/
