#include <Rcpp.h>
using namespace Rcpp;

//' het_finder function
//' 
//' Searches through the columns of a NumericMatrix and finds all the columns that 
//' contain different elements. Returns the indices of these columns. Note that 
//' this function uses 1 indexing to suit R. 
//' 
//' @param G: A NumericMatrix
//' @return A NumericVector containing the indices of the columns with different elements.
//' @examples  G <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4), het_finder(G)
//' @export
// [[Rcpp::export]]
NumericVector het_finder(NumericMatrix G) {
  int snp = G.ncol();
  int nsam = G.nrow();
  NumericVector het_indices(snp);
  std::fill(het_indices.begin(), het_indices.end(), R_NaN);
  int index = 0;
  
  //loop over all the columns of G
  for(int i=0; i<snp; i++){
    NumericVector c = G(_,i);
    bool flag = true;
    int j = 1;
    
    //loop over all elements of a column, stopping when a different element is found
    while(flag == true && j < nsam){
      if(c[0]!=c[j]){
        flag = false;
        het_indices[index] = i;
        index = index + 1;
      }
      j = j + 1;
    }
  }
  NumericVector output = het_indices[Rcpp::Range(0,index-1)] + 1;
  return output;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
