#include <Rcpp.h>
using namespace Rcpp;
#include "winsplit_length.h"
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
  //check inputs. Allow for NAs as this is not a summary statistics function.
  if(any_sug(pos<0) || any_sug(pos>1)){
    stop("Elements of pos must be between 0,1.");
  }
  if(SNP!=pos.size()){
    stop("Dimensions of pos must be the same as the number of columns in G");
  }
  if(is_sorted(pos ,true) == false){
    stop("The elements of the position vector are not in ascending order.");
  }
  
  //compute length of window
  int last_index = pos.length()-1; //C indices start at 0
  double full_len = pos[last_index]- pos[0];
  double len = full_len/n;
  
  NumericVector bounds (n+1); //initialise elements at 0.
  double start_pos = pos[0];
  //Compute the starting indices for each block
  for(int i=1;i<(n+1);i++){
    double target=start_pos+i*len;
    // Cpp starts indices at 0
    bounds[i] = first_over(pos,target) - 1; 
  }
  
  Rcout<<bounds<<std::endl;
  
  //split matrix into blocks. Each block starts at one of the boundary points. 
  //Last block ends on the final boundary points. 
  List windows(n);
  int start=bounds[0];
  int nsam=G.nrow();
  
  for(int i=0;i<n;i++){
//    Rcout<<"i is "<<i<<std::endl;
    int end=bounds[i+1] - 1;
     Rcout<<"start is"<<start<<std::endl;
     Rcout<<"end is"<<end<<std::endl;
    //if start=end. There are no SNPs in that genome window.
    if(end==bounds[i]){
      warning("No SNPs found within a window. Setting window to NULL");
      NumericMatrix m(1,1);
      std::fill(m.begin(),m.end(), NumericVector::get_na());
      windows[i]= m;
      continue;
    }
    windows[i]=G(Range(0,nsam-1),Range(start,end));
    start= bounds[i+1];
   // Rcout<<start<<std::endl;
  }
  
  return windows;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(1)
seq <-matrix(sample(0:1, size = 100, replace = TRUE), nc = 10)
pos<-runif(10, 0, 1) 
winsplit_base(seq,pos,n=4)
*/
