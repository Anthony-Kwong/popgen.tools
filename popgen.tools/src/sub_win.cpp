#include <Rcpp.h>
using namespace Rcpp;

//' sub_win function
//' 
//' Breaks down a matrix into a series of equal sized, non-overlapping windows.
//' 
//' @param G: A binary matrix of 0's and 1's. 
//' @param num_win: Number of subwindows to break G into. 
//' @return A NumericMatrix list of the windows
//' @examples sub_win(G,4)
//' @export
// [[Rcpp::export]]
List sub_win(NumericMatrix G,int num_windows) {

  int SNP=G.ncol();
  int width=SNP/num_windows;
  int nsam=G.nrow();
  
  List sub_win(num_windows);
  
  //index keeps track of which columns have been copied over
  int start=0;
  
  for(int i=0;i<num_windows;i++){
    
    if(i<num_windows-1){
      int end=start+width;
      //Rcout<<"start "<<start<<" end"<<end<<std::endl;
      sub_win[i]=G(Range(0,nsam-1), Range(start,end));
      start=end+1;
    } else {
      //Rcout<<"start "<<start<<" end"<<SNP<<std::endl;
      //remember indices start at 0
      sub_win[i]=G(Range(0,nsam-1), Range(start,SNP-1));
    }
  }

  return sub_win;
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

*/
