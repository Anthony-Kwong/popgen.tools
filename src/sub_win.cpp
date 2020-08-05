#include <Rcpp.h>
using namespace Rcpp;

//' sub_win function
//' 
//' Breaks down a matrix into a series of equal sized, non-overlapping windows. When
//' the number of columns is not divisible by the number of subwindows, the remaining
//' columns go onto the last window. 
//' 
//' @param G: A binary matrix of 0's and 1's. 
//' @param num_win: Number of subwindows to break G into. 
//' @return A List of 2 elements. 1. A list of NumericMatrix for the windows. 
//' 2. A NumericVector of column indices, indicating the boundary points for each window. 
//' We use the R setup where indices start at 1.
//' @examples sub_win(G,4)
//' @export
// [[Rcpp::export]]
List sub_win(NumericMatrix G,int num_windows) {

  int SNP=G.ncol();
  int width=SNP/num_windows; //this takes the floor because of int
  int nsam=G.nrow();
  NumericVector indices(num_windows + 1); 
  indices[0] = 0;
  
  //When there are are sub windows than SNPs, we return NULL. 
  if(width<1){
    return NULL;
  }
  
//  Rcout<<"SNP:"<<SNP<<" width:"<<width<<" nsam:"<<nsam<<std::endl;
  
  List sub_win(num_windows);
  
  //if the G.ncol() is not divisible by num_windows, we put the extra SNPs onto the last window
  
  int start=0;
  for(int i=0;i<num_windows;i++){
    
    if(i<(num_windows-1)){
      int end=start+width-1;
//      Rcout<<"start "<<start<<" end"<<end<<std::endl;
      sub_win[i]=G(Range(0,nsam-1), Range(start,end));
      indices [i+1] = end; 
      start=end+1;
    } else {
//      Rcout<<"start "<<start<<" end"<<(SNP-1)<<std::endl;
      //remember indices start at 0
      sub_win[i]=G(Range(0,nsam-1), Range(start,SNP-1));
      indices [i+1] = SNP-1;
    }
  }
  //add 1 because indices start at 1 in R
  indices = indices + 1;
  
  List output = List::create(Named("windows") = sub_win, _["bounds"] = indices);
  //List L = List::create(Named("name1") = v1 , _["name2"] = v2);
  

  return output;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(1989)
SNP=10
seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
test=sub_win(seq,3)

####

data<-readRDS("~/work/MPhil/data/toy_set.rds")
sim<-data[[213]]
seq<-sim$genomes
win_split=10
sub_win(seq,win_split)
#test1<-sum_stats(sim,win_split,100)

*/
