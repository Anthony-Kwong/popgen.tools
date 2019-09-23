##R version of theta_h for testing purposes

R_theta_h<-function(M){
  x<-colSums(M)
  N<-nrow(M)
  Si=rep(0,N-1)
  
  for(i in 1:(N-1)){
    Si[i]=count(x,i)
  }
  
  top=0
  
  for(i in 1:(N-1)){
    top=top+2*Si[i]*(i^2)
  }
  
  return (top/(N*(N-1)))
}

test_that("Fay and Wu's H computed correctly",{
  
  #count function
  
  x<-c(1,2,3,3,8,9,4,5)
  expect_equal(count(x,3),2)
  expect_equal(count(x,6),0)
  expect_equal(count(x,5),1)
  
  
  #theta_h function
  
  set.seed(2018)
  SNP=4
  seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
  expect_equal(R_theta_h(seq),theta_h(seq))

  set.seed(12344)
  SNP=50
  seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
  expect_equal(R_theta_h(seq),theta_h(seq))
  
  #fwh function
  a=14
  b=19
  expect_equal(fwh(a,b),abs(a-b))
  
})

