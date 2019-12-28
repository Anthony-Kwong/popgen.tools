#If you ever get Error in .Call(<pointer: 0x0>, N) : NULL value passed as symbol address.
#Try clearing your environment and rebuilding. 

#This function was written to test the theta_t function. 

#Count the number of pair wise differences between rows of a numeric matrix, normalised 
#by the number of pairs available. 

#input: Binary matrix G consisting of 1's and 0's. 
#output: Number of pairwise differences normalised by the number of pairs. 

test_theta_t<-function(G){
  
  pair_diff=dist(G,method="manhattan")
  total_diff=sum(pair_diff)
  nsam=nrow(G)
  total_pairs=choose(nsam,2)
  
  return(total_diff/total_pairs)
}


test_that("Tajima's D computed correctly", {

  #a1f
  N=4
  expect_equal(a1f(N),1+1/2+1/3)
  
  #a2f
  test=1+1/2^2+1/3^2
  expect_equal(a2f(N),test)
  
  #b1f
  test=(1/3)*(N+1)/(N-1)
  expect_equal(b1f(N),test)
  
  #b2f
  test=2*(N^2+N+3)/(9*N*(N-1))
  expect_equal(b2f(N),test)
  
  #c function tests
  a1=3
  a2=4
  b1=2
  b2=5
  
  #c1f
  expect_equal(c1f(b1=b1,a1=a1),b1-(1/a1))
  
  #c2f
  test=b2-(N+2)/(a1*N)+a2/(a1*a1)
  expect_equal(c2f(a1=a1,a2=a2,N=N,b2=b2),test)
  
  ##e function test
  a1=3
  a2=6
  c1=4
  c2=7
  expect_equal(e1f(a1=a1,c1=c1),c1/a1)
  expect_equal(e2f(a1=a1,c2=c2,a2=a2),c2/(a1^2+a2))
})

test_that("Tajima variance term computer correctly",{
  set.seed(1688)
  SNP=15
  seq <-matrix(sample(0:1, size =SNP*5 , replace = TRUE), nc = SNP) 
  
  nsam=nrow(seq)
  
  #Compute tajima'D coefficients
  a_1=a1f(nsam)
  a_2=a2f(nsam)
  
  b_1=b1f(nsam)
  b_2=b2f(nsam)
  
  c_1=c1f(b_1,a_1)
  c_2=c2f(a_1,a_2,b_2,nsam)
  
  e_1=e1f(c_1,a_1)
  e_2=e2f(a_1,a_2,c_2)
  
  var=e_1*nsam+e_2*nsam*(nsam-1)
  
  expect_equal(var,var_taj(seq))
})

test_that("theta_t function work correctly",{
  set.seed(2019)
  seq <-matrix(sample(0:1, size = 9, replace = TRUE), nc = 3) 
  #4 pairwise diff
  expect_equal(theta_t(seq),test_theta_t(seq))
  
  set.seed(1729)
  SNP=20
  seq <-matrix(sample(0:1, size = SNP*8, replace = TRUE), nc = SNP) 
  expect_equal(theta_t(seq),test_theta_t(seq))
  
  
})

test_that("Tajima D computed correctly",{
  set.seed(1762)
  
  nums=rnorm(3,mean=5,sd=3)
  
  t=nums[1]
  w=nums[2]
  v=nums[3]
  
  test=(t-w)*v^(-0.5)
  
  expect_equal(taj_D(t,w,v),test)
  
})


