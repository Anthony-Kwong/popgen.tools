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
  //Rcout << sites <<std::endl;
  int nsam = G.nrow();
  // loop across the columns of G
  for(int c = 0; c < sites; c++){
    //Rcout << "c " << c << std::endl; 
    //copy the c^th column
    NumericVector col = G(_,c);
    //find all elements in column which are NA
    NumericVector na_pos = which(f_na(col));
    //the which function from R uses 1 indexing. Converting to 0 indexing for cpp.
    na_pos = na_pos - 1; 
    int num_na = na_pos.length();
    
    //if we run into a column with only NAs, convert that column to a vector of 0's.
    if(num_na==nsam){
      output(_,c) = NumericVector(nsam);
      continue;
    }
    
    // obtain all the non-NA elements in the c^th column
    col = discard(col, f_na);
   // Rcout << col << std::endl;
    
    //loop over all NAs and impute them
    for(int r=0; r<num_na; r++){
      
      //note r1 is the first NA, not the first row
      //Rcout <<"r "<< r << std::endl;
      
      //random sample from non-NA elements
      NumericVector imp = runif(1,Named("min")=0,_["max"]=col.length());
      //round down to nearest integer to get a valid index
      int imp_value = imp[0];
      // Rcout <<"imp "<< imp <<std::endl;
      // Rcout << "imp_value "<< imp_value <<std::endl;
      // Rcout << "change " << col[imp_index]<<std::endl;
      // Rcout << "pos_na"<< na_pos << std::endl;
      // Rcout<<"num_na"<< num_na <<std::endl; 
      
      //change the NA element
      output(na_pos[r],c) = col[imp_value];
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
