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

//' vec_split function
//' 
//' Breaks a vector x into n equal sized blocks. If the dimension of x is not divisible by n, we put the remaining elements into the last block. 
//' 
//' @param x: vector to split
//' @param n: number of blocks to split a vector into. 
//' @return a list of vectors
//' @examples vec_split(seq(10),3)
//' @export
// [[Rcpp::export]]
Rcpp::List vec_split(NumericVector x, int n) {
  int dim=x.length();
  double len=1.0*dim/n;
  int width=len;
  
  //len=round(len);
  List vec_blocks(n);
  
  //add start and end indices
  int start=0;
  int end=width-1;
    
  for(int i=0;i<n;i++){
    if(i<n-1){
      //standard filling of blocks
      vec_blocks[i]=x[Rcpp::Range(start,end)];
      start=end+1;
      end=end+width;
    } else {
      //putting remainder elements into the last block
//      Rcout<<"trigger"<<std::endl;
//      Rcout<<"last start "<<x[start]<<std::endl;
//      Rcout<<"last one "<<x[dim-1]<<std::endl;
      vec_blocks[i]=x[Rcpp::Range(start,dim-1)];
    }
  }
  
  
  //for debugging purposes. 
  
  // Rcout<<"width is "<<width<<std::endl;
  // Rcout<<"list has "<< n<<" elements"<<std::endl;
  // Rcout<<"dim of vector is "<<dim<<std::endl;
  // 
  // for(int i=0;i<n;i++){
  //   //Rcout<<vec_blocks[i]<<std::end;
  //   Rf_PrintValue(vec_blocks[i]);
  // }
  
  
  return vec_blocks;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
x<-seq(20)
vec_split(x,5)

x<-seq(7)
vec_split(x,2)
*/
