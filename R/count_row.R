count_row<-function(A,x){
  count=0
  for(i in 1:(nrow(seq))){
    if(identical(A[i,],as.integer(x))==TRUE){
      count=count+1
    }
  }
  return (count)
}