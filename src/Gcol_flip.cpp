#include <Rcpp.h>
#include "Gcol_flip.h"
using namespace Rcpp;

//' Gcol_flip function
//' 
//' Takes binary genome matrix and a column index. Returns the number of pairwise
//' differences that would occur if the individual elements of the column were 
//' flipped. 
//' 
//' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
//' @param col: A column index. 0 indexing is used. 
//' @return A NumericVector where element i is the number of pairwise differences in G if 
//' element i in the selected column was flipped. 
//' @export
// [[Rcpp::export]]
NumericVector Gcol_flip(NumericMatrix G, int col) {
  int nsam = G.nrow();
  int snp = G.ncol();
  // initialise vector to carry the number of pairwise differences
  NumericVector pd(nsam);
  
  //check input
  if(col > (snp - 1)){
    Rcpp::stop("input col index exceeds number of columns.");
  }
  if(col < 0){
    Rcpp::stop("input col index is smaller than 0.");
  }
  
  //loop over each column
  for(int i=0; i < nsam; i++){
    //make deep copy of matrix over to modify. 
    Rcpp::NumericMatrix K(Rcpp::clone(G));
    int el = K(i,col);
    
    //flip the element
    if(el==1){
      //Rcout<<"flip 1"<<std::endl;
      K(i,col) = 0;
    } else if (el==0){
      //Rcout<<"flip 0"<<std::endl;
      K(i,col) = 1;
    } else {
      Rcpp::stop("Error in Gcol_flip. The elements of the input matrix must be either 0 or 1");
    }
    //Rcout<<K<<std::endl;
    //compute the number of pairwise differences
    pd(i) = pairwise_diff(K);
  }
  
  return pd;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
