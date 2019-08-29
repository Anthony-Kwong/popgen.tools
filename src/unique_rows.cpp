#include <Rcpp.h>
using namespace Rcpp;
#include "vec_equal.h"

//including the header allows us to use that function on this script

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' fill_row function
//' takes a NumericVector and adds it as a new row on the bottom of a matrix
//' @param A: A NumericMatrix.
//' @param x: A NumericVector. 
// [[Rcpp::export]]
NumericMatrix fill_row(NumericMatrix A, NumericVector x){
  int colA=A.ncol();
  int size_x=x.size();
  
  if(colA!=size_x){
    try{
      throw 1;
      }
    catch(int e){
      Rcout<<"Error exception "<<e<<std::endl;
      Rcout<<"Dimensions of vector must match the number of columns in NumericMatrix"<<std::endl;
      //throw an null matrix
      NumericMatrix Z;
      return Z;
    }
  }
    
//    Rcout<<t<<std::endl;
    NumericMatrix B (A.nrow()+1, A.ncol());
    for(int i=0;i<A.nrow();i++){
      B(i,_)=A(i,_);
    }
    //insert vector as the last row of matrix
    B(A.nrow(),_)=x;
  return B;
}

//make a check function. takes a vector and checks if it is a row in a Matrix

//' present_row function
//' 
//' checks if a vector is present as a row in the input matrix
//' 
//' @param A: NumericMatrix
//' @param x: NumericVector
//' @return boolean: true if the vector is present as a row in the matrix
//' @export
// [[Rcpp::export]]
bool present_row(NumericMatrix A, NumericVector x){
  int ncols=A.ncol();
  
  if(ncols!=x.size()){
    //error message
    try{
      throw 1;
    }
    catch(int e){
      Rcout<<"present_row error. Exception "<<e<<std::endl;
      Rcout<<"Error: Dimensions of vector x must match the number of columns in matrix A"<<std::endl;
      return FALSE;
    }
  }
  
//  Rcout<<A.nrow()<<std::endl;
  int num_rows=A.nrow();
  for(int i=0;i<num_rows;i++){
    NumericVector row=A(i,_);
//    Rcout<<i<<std::endl;
//    Rcout<<row<<std::endl;
    if(vec_equal(row,x)==TRUE){
//      Rcout<<"We found one folks!"<<std::endl;
      return TRUE;
    }
  }
  return FALSE;
}

//modify to return a list, keep track of haplotype frequencies

//' unique_rows function
//' 
//' Takes a NumericMatrix and returns all the unique rows as a NumericMatrix
//' 
//' @param A: A general matrix of real values.
//' @return A NumericMatrix whose rows are all the unique rows of input matrix A.
//' @examples unique_rows(A)
//' @export
// [[Rcpp::export]]
NumericMatrix unique_rows(NumericMatrix A) {
  //intialise matrix with dim nrow*nncol
  NumericMatrix B(A.nrow(),A.ncol());
  std::fill( B.begin(), B.end(), NumericVector::get_na());
//  Rcout<<B<<std::endl;
  
  
  //index keeps track of how many rows we have copied over. 
  int index=0;
  
  //change design, keep indices of unique rows
  
  //copy first row of A into B
  NumericVector row=A(0,_);
  B(0,_)=row;
  
  //loop across all rows
  for(int i=1; i<A.nrow(); i++){
    row=A(i,_);
    if(present_row(B,row)==FALSE){
      B(i,_)=row;
      // Rcout<<"hit"<<std::endl;
      // Rcout<<row<<std::endl;
      // Rcout<<B<<std::endl;
      index+=1; 
    }
  }
//  Rcout<<B<<std::endl;
  NumericMatrix C= B(Range(0,index),Range(0,B.ncol()-1));
//  Rcout<<C<<std::endl;
  
  
  return C;
}







// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(1709)
seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4) 
x<-c(1,0,1,1)
seq<-fill_row(seq,x)

unique_rows(seq)
*/
