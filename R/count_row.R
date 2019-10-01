#This function was written for testing purposes only. It is not part of the main package.

#count_row function: Count the number of times vector x appears as a row in matrix A. 

#inputs: A: an m by n matrix, x: an n dimensional vector
#output: scalar output of the number of times x appears as a row in A. 
count_row<-function(A,x){
  count=0
  for(i in 1:(nrow(seq))){
    if(identical(A[i,],as.integer(x))==TRUE){
      count=count+1
    }
  }
  return (count)
}