test_that("row_replicate works",{
  set.seed(32)
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)
  
  N=5
  index=2
  output<-row_replicate(seq,row=index,n=N)
  
  for(i in 1:N){
    seq=rbind(seq,seq[index,]) 
  }
  
  test=identical(seq,output)
  expect_equal(test,TRUE)
  
})