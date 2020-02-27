#include <Rcpp.h>
using namespace Rcpp;
#include "winsplit_length.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' winsplit_base function
//' 
//' Split a genome matrix into equally sized windows as determined by
//' the number of bases. The windows can differ in the number of columns
//' because SNPs are not always uniform across a given region. 
//' 
//' @param G: A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. Each column is a SNP. Each row is a sampled individual. 
//' @param pos: A NumericVector consisting of values between 0,1. The ith element is the position of the ith SNP in G.
//' @param n: integer for the number of desired windows in the output list. 
//' @return A NumericMatrix list of the windows
//' @examples seq <-matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
//' pos<-runif(0,1,n=5) %>% sort()
//' winsplit_length(seq,pos,5)
//' @export
// [[Rcpp::export]]
List winsplit_base(NumericMatrix G, NumericVector pos, int n) {
  
  int SNP=G.ncol();
  //check inputs
  if(any_sug(pos<0) || any_sug(pos>1)){
    stop("Elements of pos must be between 0,1.");
  }
  if(SNP!=pos.size()){
    stop("Dimensions of pos must be the same as the number of columns in G");
  }
  
  double len = 1.0/n;
  NumericVector start_indices(n+1);
  
  for(int i=1;i<(n+1);i++){
//    Rcout<<i<<std::endl;
    start_indices[i] = find_index(pos,i*len);
  }
  
//  Rcout<<start_indices<<std::endl;
  
  List windows(n);
  int start=start_indices[0];
  int nsam=G.nrow();
  
  for(int i=0;i<n;i++){
    int end=start_indices[i+1];
    Rcout<<"start is"<<start<<std::endl;
    Rcout<<"end is"<<end<<std::endl;
    //if start=end. There are no SNPs in that genome window.
    if(end==start_indices[i]){
      warning("No SNPs found within a window. Setting window to NULL");
      windows[i]=NULL;
      break;
    }
    windows[i]=G(Range(0,nsam-1),Range(start,end));
    start=end+1;
   // Rcout<<start<<std::endl;
  }
  
  return windows;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
seq <-matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
pos<-seq(0,1,by=0.2)
winsplit_base(seq,pos,n=4)
*/
