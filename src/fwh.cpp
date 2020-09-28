#include <Rcpp.h>
using namespace Rcpp;
#include "matrix_check.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' count function
//' 
//' counts the number of times a value appears in an IntegerVector
//' @param vec: a NumericVector
//' @param target: input integer value to find
//' @return scalar value of the number of times the target appears in vec
//' @examples count(c(1,2,3,4),3)
//' @export
//[[Rcpp::export]]
int count(NumericVector vec,int target){
  int counter=0;
  for(int i=0;i<vec.size();i++){
    if(vec[i]==target){
      counter+=1;
    }
  }
  return counter;
}



//' theta_h function
//' 
//' Computes theta_h for a genome matrix
//' 
//' @param G: Binary genome matrix with 0's and 1's. Each column is a SNP, each row is an individual.
//' @return scalar value of theta_h for that sampled population. 
//' @examples theta_h(G)
//' @export
// [[Rcpp::export]]
double theta_h(NumericMatrix G){
  
  //check input
  if(have_na(G)){
    Rcpp::warning("Input matrix is NA. Returning NA.");
    return (R_NaN);
  }
  
  //removing feature because we still want to compute a value for aDNA sims with few columns
  
  // if(G.ncol()<5){
  //   Rcpp::warning("Less than 5 SNPs in genome matrix. Returning NA.");
  //   return (R_NaN);
  // }
  
  //compute the frequency of each SNP/mutation
  NumericVector column_sum = colSums(G);
  
  //take the max frequency to make sum efficient
  int max_derived = max(column_sum);

  //S_i is the number of mutations with frequency i among the sampled sequences. 
  //Preallocate memory for max+1 since cpp indices start at zero. The S_i[0]=0 and won't get used in the sum. 
  NumericVector S_i(max_derived + 1);
  S_i[0]=0;

  for(int i = 1;i <= max_derived; i++){
    S_i[i] = count(column_sum,i);
  }
 
 //compute the sum term. 
  
  double sum_term=0;
  
  for(int i = 1;i <= max_derived; i++){
   double temp = S_i[i]*i*i;
   sum_term = sum_term+temp;
  }
  
  int N = G.nrow();
  double H = sum_term/(N*(N-1)/2);
  return H;
}

//' fwh function
//' 
//' Computes Fay and Wu's H for a genome matrix. Strength of H indicates magnitude of selective sweep. H=0 indicates there is no evidence of deviation from neutrality. See doi: 10.1534/genetics.106.061432.
//' 
//' @param t_t: theta_t for genome matrix G. Use theta_t(). Also called theta_pi in literature, 
//' @param t_h: theta_h for genome matrix G. Use theta_h().
//' @return scalar value of Fay and Wu's H for that sampled population.  
//' @examples fwh(theta_t(G),theta_h(G))
//' @export
// [[Rcpp::export]]
double fwh(double t_t,double t_h) {
  double stat=t_t-t_h;
  return stat;
}

//H is the difference between theta_pi (average het) and theta_h

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(800)
seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)
theta_h(seq)

count(c(1,2,3,3),4)

*/
