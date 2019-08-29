test_that("fill_row function works",{
  set.seed(1492)
  seq <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 4) 
  x<-c(1,0,1,1)
  test1<-fill_row(seq,x)
  test2<-rbind(seq,x,deparse.level = 0)
  expect_equal(test1,test2)
})

test_that("present_row function works",{
  set.seed(1492)
  seq <-matrix(sample(0:1, size = 80, replace = TRUE), nc = 8)
  expect_equal(present_row(seq,seq[4,]),TRUE)
  
  set.seed(1866)
  seq <-matrix(sample(0:1, size = 96, replace = TRUE), nc = 8)
  expect_equal(present_row(seq,seq[9,]),TRUE)
  
  set.seed(1865)
  seq <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 8)
  expect_equal(present_row(seq,c(0,0,0,0,0,0,0,1)),FALSE)
})