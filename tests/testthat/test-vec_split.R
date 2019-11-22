test_that("vec_split computed correctly",{
  
  #no remaining elements case
  x<-seq(12)
  n<-3
  output<-vec_split(x,n)
  actual<-list(seq(1,4),seq(5,8),seq(9,12))
  expect_equal(output,actual)
  
  #remaining elements case. Want extra elements to go to the last block. 
  
  y<-seq(15)
  n<-4
  output<-vec_split(y,n)
  actual<-list(seq(1,3),seq(4,6),seq(7,9),seq(10,15))
  expect_equal(output,actual)
})