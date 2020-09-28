##R version of theta_h for testing purposes

R_theta_h<-function(M){
  
  if(ncol(M)<5){
    return (NaN)
  }
  x<-colSums(M)
  N<-nrow(M)
  Si=rep(0,N-1)
  
  for(i in 1:(N-1)){
    Si[i]=count(x,i)
  }
  
  sum_term=0
  
  for(i in 1:(N-1)){
    sum_term=sum_term+Si[i]*(i^2)
  }
  
  return (sum_term/choose(N,2))
}

test_that("Fay and Wu's H computed correctly",{
  
  #count function
  
  x<-c(1,2,3,3,8,9,4,5)
  expect_equal(count(x,3),2)
  expect_equal(count(x,6),0)
  expect_equal(count(x,5),1)
  
  
  #theta_h function
  
  #when inputing sample matrices, make sure that you don't get entire columns of 1's. This does not happen in actual genome simulations.
  
  set.seed(2018)
  SNP=4
  seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
  suppressWarnings(
    expect_equal(R_theta_h(seq),theta_h(seq))
  )

  set.seed(123)
  SNP=50
  seq <-matrix(sample(0:1, size = SNP*5, replace = TRUE), nc = SNP)
  expect_equal(R_theta_h(seq),theta_h(seq))
  
  set.seed(1234)
  SNP=500
  seq <-matrix(sample(0:1, size = SNP*20, replace = TRUE), nc = SNP)
  #Just checking that we don't have any columns with just 1's. 
  #colSums(seq) %>% max()
  expect_equal(R_theta_h(seq),theta_h(seq))
  
  # #test NA input returns NA. 
  # null_win = as.numeric(NA) %>% matrix()
  # check = suppressWarnings(  theta_h(null_win)  )
  # expect_equal(check, NaN)
  # 
  # #fwh function
  # a=14
  # b=19
  # expect_equal(fwh(a,b),a-b)
  # expect_equal(fwh(NaN,NaN), NaN)
  
})

