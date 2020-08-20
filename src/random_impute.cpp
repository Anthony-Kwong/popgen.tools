#include <Rcpp.h>
using namespace Rcpp;

//' random_impute function
//' 
//' Imputes NAs in a genome matrix. For each NA, we randomly sample from the non-NA
//' elements of that column. Each element is equally likely to be selected. NA becomes
//' the randomly selected element. 
//' 
//' @param G: A NumericMatrix with NAs. 
//' @param seed: A numeric random seed.
//' @return A NumericMatrix with imputed values in place of the NAs. 
//' @examples G <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)
//' G[1,1] = NA
//' random_impute(G)
//' @export
// [[Rcpp::export]]
NumericMatrix random_impute(NumericMatrix G) {
  Environment pkg = Environment::namespace_env("purrr");
  Function discard = pkg["discard"];
  Function f_na("is.na");
  Function which("which");
  Function runif("runif");
  
  NumericMatrix output = G;
  int sites = G.ncol();
  int nsam = G.nrow();
  for(int c=0; c < sites; c++){
    NumericVector col = G(_,c);
    NumericVector na_pos = which(f_na(col));
    na_pos = na_pos - 1; //the which function uses 1 indexing. Converting to 0 indexing for cpp.
    int num_na = na_pos.length();
    col = discard(col, f_na);
    for(int r=0; r<num_na; r++){
      //Rcout << r << std::endl;
      //random sample from col
      NumericVector imp = runif(1,Named("min")=0,_["max"]=col.length());
      int imp_index = imp[0];
      Rcout <<"imp "<< imp <<std::endl;
      Rcout << "imp_index "<< imp_index <<std::endl;
      Rcout << "change "<< col[imp_index]<<std::endl;
      // Rcout << "pos_na"<< na_pos << std::endl;
      // Rcout<<"num_na"<< num_na <<std::endl; 
      output(na_pos[r],c) = col[imp_index];
    }
  }
  return output;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
