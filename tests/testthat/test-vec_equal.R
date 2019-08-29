test_that("vec_equal function works",{
  
  #small case
  x<-c(1,0,0)
  y<-c(1,0,0)
  expect_equal(vec_equal(x,y),TRUE)
  
  y<-c(1,1,1)
  expect_equal(vec_equal(x,y),FALSE)
  
  #large case
  set.seed(1066)
  seq1<-matrix(sample(0:1, size = 50, replace = TRUE), nc = 50) 
  
  set.seed(1066)
  seq2<-matrix(sample(0:1, size = 50, replace = TRUE), nc = 50) 
  
  expect_equal(vec_equal(seq1,seq2),TRUE)
})