#test function for winsplit_base
#
#Breaks a genome matrix into equal chunks based on genome length
#
#param: NumericMatrix G, binary genome matrix
#param: NumericVector pos, sorted vector with values between 0,1
#param: integer n. number of chunks to make.
#return: A list of matrices.

R_winsplit_base<-function(G, pos, n){
  length = 1/n
  wins = list()
  start_indices=rep(NA,n+1)
  
  start_indices[1]=1
  for(i in 2:(n+1)){
    dist=pos-(i-1)*length
    dist=abs(dist)
    index=which.min(dist)
    start_indices[i]=index
  }
  
  start=1
  for(j in 2:(n+1)){
    end=start_indices[j]
    wins[[j-1]]=G[,start:end]
    start=end+1
  }
  
  return(wins)
}


test_that("winsplit_base works",{
  #manual check
  set.seed(20)
  SNP=10
  seq <-matrix(sample(0:1, size = SNP*6, replace = TRUE), nc = SNP)
  pos <- runif(0,1,n=SNP) %>% sort()
  C_says <- winsplit_base(seq,pos,n=3)
  R_says <- R_winsplit_base(seq,pos,n=3)
  
  ans=list()
  ans[[1]] = seq[,1:5]
  ans[[2]] = seq[,6:7]
  ans[[3]] = seq[,8:10]
  
  #R version
  for(i in 1:3){
    check=all.equal(R_says[[i]],ans[[i]])
    expect_equal(check,T)
  }
  
  #C version
  for(i in 1:3){
    check=all.equal(C_says[[i]],ans[[i]])
    expect_equal(check,T)
  }
  
  #autocheck
  
  set.seed(21)
  SNP=250
  nwins=10
  seq <-matrix(sample(0:1, size = SNP*50, replace = TRUE), nc = SNP)
  pos <- runif(0,1,n=SNP) %>% sort()
  C_says <- winsplit_base(seq,pos,n=nwins)
  R_says <- R_winsplit_base(seq,pos,n=nwins)
  
  for(i in 1:nwins){
    check=all.equal(C_says[[i]],R_says[[i]])
    expect_equal(check,T)
  }
})
