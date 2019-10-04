#' row_replicate function
#' 
#' Replicates a selected row in a matrix n times and appends it to the bottom. 
#'
#' @param M: A numeric matrix
#' @param row: the index of the row  to replicate
#' @param n: number of times to replicate row
#'
#' @return A numeric matrix with selected row appended at the bottom
#' @export
#'
#' @examples row_replicate(M,row=2,n=3)
row_replicate<-function(M,row,n){
  r=M[row,]
  M2=M
  for(i in 1:n){
    M2=rbind(M2,r,deparse.level = 0)
  }
  
  return (M2)
}
