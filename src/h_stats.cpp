#include <Rcpp.h>
using namespace Rcpp;
#include "unique_rows.h"
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

//remember param is singular for the documentation

//' h_stats function
//' 
//' Computes the h1,h2,h12,h123 statistics. For more information on how these are computed see https://doi.org/10.1371/journal.pgen.1005004.
//' 
//' @param G: A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. Each column is a SNP. Each row is a sampled individual. 
//' @return A numeric vector of h stats. (h1,h2,h12,h123)
//' @examples G = matrix(sample(0:1, size =5*5 , replace = TRUE), nc = 5), h_stats(G)
//' @export
// [[Rcpp::export]]
NumericVector h_stats(NumericMatrix G) {
  //check input
  if(have_na(G)){
    Rcpp::warning("Input matrix is NA. Returning NA.");
    NumericVector h_stats (4);
    std::fill(h_stats.begin(), h_stats.end(), NumericVector::get_na());
    return (h_stats);
  }
  
  //removing feature because we still want to compute a value for aDNA sims with few columns
  
  // if(G.ncol()<5){
  //   Rcpp::warning("Less than 5 SNPs in genome matrix. Returning NA.");
  //   NumericVector h_stats (4);
  //   std::fill(h_stats.begin(), h_stats.end(), NumericVector::get_na());
  //   return (h_stats);
  // }
  
  //compute the haplotype frequencies
  NumericVector freq = row_freq(G);
  
  //sort the frequencies in descending order
  NumericVector p = vec_sort(freq);
  
  int num_haplotypes=p.size();

  //compute h1. Sum of squares of all haplotype frequencies. 
  double h1=0;
  for(int i=0; i<num_haplotypes; i++){
    h1=h1+p[i]*p[i];
  }
   
  //the h12 stat combines the frequencies of the top 2 most common haplotypes into one haplotype. We use a shortcut to compute it. 
  double h12=h1+2*p[0]*p[1];

  //h3 stat is like h2 but combines the top 3. We use algebraic shortcut.
  double h123=h12+2*p[0]*p[2]+2*p[1]*p[2];

  //h2. This is like h1 but removes the first haplotype.
  double h2=h1-p[0]*p[0];

  //Putting everything together
  NumericVector h_stats=NumericVector::create(h1,h2,h12,h123);

  return h_stats;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(800)
seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)
seq<-rbind(seq,seq[2,])
seq<-rbind(seq,seq[4,])
seq<-rbind(seq,seq[3,])
h_stats(seq)
*/
