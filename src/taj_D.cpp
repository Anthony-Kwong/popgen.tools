#include <Rcpp.h>
using namespace Rcpp;


//////////////////////////////////////////////////////////////////////////////

//functions to compute the coefficients for Tajima's D. These require the number of
// sequences (N), number of segregating sites (nseg) and average pairwise differences
// (pi) to be known.


//Call Tajima last with header
// [[Rcpp::export]]
double a1f(int N){
  double a1;
  a1=0.0;
  int i;
  for(i=1; i<=N-1; i++){
    a1=a1+(1.0/i);
    // std::cout << "loop" << std::endl;
    
    //    std::cout << a1 << std::endl;
  }
  return (a1);
}

// [[Rcpp::export]]
double a2f(int N){
  double a2;
  a2=0.0;
  int i;
  for(i=1; i<=N-1; i++){
    a2+=1.0/(i*i);
  }
  return (a2);
}

// [[Rcpp::export]]
double b1f(double N){
  double b1;
  b1=(N+1)/(3*(N-1));
  Rcout<<b1<<std::endl;
  return(b1);
}

// [[Rcpp::export]]
double b2f(double N){
  double b2;
  double top=2*(N*N+N+3);
  double bot=9*N*(N-1);
  b2=top/bot;
  return(b2);
}

// [[Rcpp::export]]
double c1f(double b1, double a1){
  double c1;
  c1=b1-(1/a1);
  return(c1);
}

// [[Rcpp::export]]
double c2f(double a1, double a2, double b2, int N){
  double c2;
  c2=b2-(N+2)/(a1*N)+a2/(a1*a1);
  return (c2);
}

// [[Rcpp::export]]
double e1f(double c1, double a1){
  double e1;
  e1=c1/a1;
  return(e1);
}

// [[Rcpp::export]]
double e2f(double a1, double a2, double c2){
  double e2;
  e2=c2/((a1*a1)+a2);
  return(e2);
}

////////////////////////////////////////////////

//Theta_t: Compute pairwise differences. In popgen, theta_t is also called pi.
//Input: G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.

//' theta_t function
//' Computes theta_t, the number of pairwise differences normalised by the number of pairs. 
//' @param G: G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
//' @return scalar value of theta_t
// [[Rcpp::export]]
double theta_t(NumericMatrix G){
  int num_sam=G.nrow();
  int num_seg=G.ncol();
  
  //add in a array that keeps track of differences in each pair,
  //length of array is the number of possible pairs.
  int total_pairs=num_sam*(num_sam-1)/2;
  NumericVector diff(total_pairs);
  
  
  //k is the index for the the pairs.
  int k=0;
  
  //loop through individuals. Remember indices start at 0 in cpp
  
  for (int n=0; n<=num_sam-1; n++){
    
    //loop through each pair
    for (int i=n+1; i<=num_sam-1;i++){
      
      //loop through each SNP
      for (int j=0; j<=num_seg-1; j++){
        
        if(G(n,j)!=G(i,j)){
          //index -1 because indices start at 0
          diff(k)+=1;
        }
      }
      //      Rcout<<"Pair"<<k<<" diff"<<diff(k)<<std::endl;
      k+=1;
    }
  }
  //we need total to be a double. C++ prioritizes the nominator in division.
  double total=sum(diff);
//  Rcout<<diff<<std::endl;
//  Rcout<<total<<" "<<total_pairs<<std::endl;
  double D=total/total_pairs;
  return (D);
}

//'theta_w function
//' Computes the average number of pairwise differences. 
//' @param G: Binary genome matrix of 0's and 1's. Each column is a SNP, each row is an individual.
//' @return scalar value opf theta_w
// [[Rcpp::export]]
double theta_w(NumericMatrix G){
  int N=G.nrow();
  int num_seg=G.ncol();
  double theta_w;
  theta_w=num_seg/a1f(N);
  return (theta_w);
}

//remember this header has to be immediately above the Rcpp::export. 
//otherwise documentation won't be generated. 

//' taj_D function
//' @param G: A binary matrix of 0's and 1's. Each column is a SNP and each row is a sampled individual.
//' @return a scalar value of tajima's D for the sampled population. 
//' @examples taj_D(G)
//' @export
// [[Rcpp::export]]
double taj_D(NumericMatrix G){
  int nsam=G.nrow();
  
  //compute Tajima_D coefficients
  double a_1=a1f(nsam);
  double a_2=a2f(nsam);
  
  double b_1=b1f(nsam);
  double b_2=b2f(nsam);
  
  //double c1f(double b1, double a1)
  double c_1=c1f(b_1,a_1);
  //double c2f(double a1, double a2, double b2, int N)
  double c_2=c2f(a_1,a_2,b_2,nsam);
  
  //double e1f(double c1, double a1)
  //double e2f(double a1, double a2, double c2)
  double e_1=e1f(c_1,a_1);
  double e_2=e2f(a_1,a_2,c_2);
  
  double var=e_1*nsam+e_2*nsam*(nsam-1);
  
  double top=theta_t(G)-theta_w(G);
  return(top/pow(var,0.5));
}


//var(d)=e1*num_seg+e2*num_seg(num_seg-1)





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
b1f(4)
*/
