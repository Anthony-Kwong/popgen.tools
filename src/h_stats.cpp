#include <Rcpp.h>
using namespace Rcpp;
#include "vec_equal.h"
#include "unique_rows.h"

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
//' Computes the h1,h2,h12 statistics. For more information see https://doi.org/10.1371/journal.pgen.1005004.
//' 
//' @param G: A binary genome matrix consisting of 1's and 0's. Each column is a SNP. Each row is a sampled individual. 
//' @return A numeric vector of h stats. (h1,h2,h12)
//' @examples h_stats(G)
//' @export
// [[Rcpp::export]]
int h_stats(NumericMatrix G) {
  NumericVector haplo_counts=unique_rows(G);
  int num_haplotypes=haplo_counts.size();
  int nsam=G.nrow();
  Rcout<<haplo_counts<<std::endl;
  //Rcout<<p<<" Before"<<std::endl;
  
  double p [num_haplotypes];
  //Rcout<<p<<std::endl;
  
  //convert counts to frequencies. I used NumericVector before because of functional support
  //and ease of debugging. 
  for(int i=0; i<num_haplotypes;i++){
    p[i]=haplo_counts[i]*(1.0/nsam);
    //Rcout<<p[i]<<std::endl;
  }
  
  // //find the frequencies of 1st, 2nd, 3rd most common haplotypes
  // //if we ever need more frequencies use std::sort
  // int index=which_max(haplo_counts);
  // double p_1=haplo_counts[index];
  // haplo_counts[index]=NumericVector::get_na();
  // 
  // int index=which_max(haplo_counts);
  // double p_2=haplo_counts[index];
  // haplo_counts[index]=NumericVector::get_na();
  
  
  
  // Rcout<<p_1<<std::endl;
  // Rcout<<haplo_counts<<std::endl;

  //p=p*(1/G.nrow());
  //Rcout<<p<<" After"<<std::endl;
  //Rcout<<p.size()<<std::endl;
  //Rcout<<G.nrow()<<std::endl;
  
  //compute h1
  double h1=0;
   for(int i=0; i<num_haplotypes; i++){
     h1=h1+p[i]*p[i];
   }
   //Rcout<<h1<<std::endl;
   
   
   
   //which_max(x)
   
   //h12 need biggest 2 frequences
   
   //h123 need biggest 3 frequencies
  
  
  return 0;
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
