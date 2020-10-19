#include <Rcpp.h>
#include "Gcol_denoise.h"
using namespace Rcpp;

//' Gcol_denoise function
//' 
//' Takes a matrix and a column index. Flips the element in that column which leads
//' to the greatest reduction in the number of pairwise differences. When there is 
//' a tie, the first element among the tied elements is changed. 
//' 
//' @param: G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
//' @param: col: A column index. 0 indexing is used. 
//' @return: A NumericMatrix like G, with one element flipped. 
//' @examples:   G = matrix(sample(0:1, size =15 , replace = TRUE), nc = 5)
//' Gcol_denoise(G,3)
//' @export
// [[Rcpp::export]]
NumericMatrix Gcol_denoise(NumericMatrix G, int col) {
  
  int snp = G.ncol();

  //check input
  if(col > (snp - 1)){
    Rcpp::stop("input col index exceeds number of columns.");
  }
  if(col < 0){
    Rcpp::stop("input col index is smaller than 0.");
  }
  
  //find which element to flip
  NumericVector x = Gcol_flip(G,col);
  
  //if all flips produce the same number of pwd, return original matrix
  int y = unique(x).length();
  if(y < 2){
    return G;
  }
  
  //which_min returns the index of the first smallest element in a vector
  int flip_index = which_min(x);
  //Rcout << flip_index << std::endl;
  
  //flip element in the matrix
  int el = G(flip_index,col);
  if(el == 0){
    G(flip_index,col) = 1;
  } else if (el == 1) {
    G(flip_index,col) = 0;
  } else {
    Rcpp::stop("Error in Gcol_flip. The elements of the input matrix must be either 0 or 1");
  }
  
  return G;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
