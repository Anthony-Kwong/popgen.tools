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
NumericVector h_stats(NumericMatrix G) {
  NumericVector haplo_counts=unique_rows(G);
  int num_haplotypes=haplo_counts.size();
  int nsam=G.nrow();
  //Rcout<<haplo_counts<<std::endl;
  //Rcout<<p<<" Before"<<std::endl;
  //
  
  //double vector to store the frequencies of each haplotype
  double p [num_haplotypes];
  //Rcout<<p<<std::endl;
  
  //flag for change. unique_rows should give frequencies already.
  
  //convert counts to frequencies. I used NumericVector before because of functional support
  //and ease of debugging. 
  for(int i=0; i<num_haplotypes;i++){
    p[i]=haplo_counts[i]*(1.0/nsam);
    //Rcout<<p[i]<<std::endl;
  }
  
  //compute h1
  double h1=0;
  for(int i=0; i<num_haplotypes; i++){
    h1=h1+p[i]*p[i];
  }
  //Rcout<<"h1 "<<h1<<std::endl;
  
  //find the frequencies of 1st, 2nd, 3rd most common haplotypes
  // NumericVector top_hap=three_top(haplo_counts);
  // double top_freqs[3];
  // for(int i=0;i<3;i++){
  //   top_freqs[i]=top_hap[i]*(1.0/nsam);
  // }
  
  //Rcout<<top_hap<<std::endl;
  
  //testing if frequencies were calculated correctly
 // for(int i=0;i<3;i++){
    //Rcout<<top_freqs[i]<<std::endl;
 // }
   
   //the h2 stat combines the frequencies of the top 2 most common haplotypes into one haplotype.
   // double h2=h1+2*top_freqs[0]*top_freqs[1];
   // //Rcout<<"h2 "<<h2<<std::endl;
   // 
   // //h3 stat is like h2 but takes top 3 
   // double h3=h2+2*top_freqs[0]*top_freqs[2]+2*top_freqs[1]*top_freqs[2];
   // //Rcout<<"h3 "<<h3<<std::endl;
   // 
   // NumericVector h_stats=NumericVector::create(h1,h2,h3);
   //Rcout<<h_stats<<std::endl;
  
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
