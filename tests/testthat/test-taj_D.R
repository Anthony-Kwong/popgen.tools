#pacman::p_load(testthat,popgen.tools)

#If you ever get Error in .Call(<pointer: 0x0>, N) : NULL value passed as symbol address.
#Try clearing your environment and rebuilding. 


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


