library("popgen.tools")

test_that("fill_row function works",{
  set.seed(1492)
  seq <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 4) 
  x<-c(1,0,1,1)
  test1<-fill_row(seq,x)
  test2<-rbind(seq,x,deparse.level = 0)
  expect_equal(test1,test2)
})

# count_row<-function(A,x){
#   count=0
#   for(i in 1:(nrow(seq))){
#     if(identical(A[i,],as.integer(x))==TRUE){
#       count=count+1
#     }
#   }
#   return (count)
# }

# test_that("row count test function works",{
#   set.seed(1492)
#   seq <-matrix(sample(0:1, size = 20, replace = TRUE), nc = 4)
#   x<-c(1,0,1,0)
#   y<-c(1,1,0,1)
#   z<-c(1,1,1,1)
#    test1<-count_row(seq,x)
# })


test_that("present_row function works",{
  #eventually fix this so we compare with R function that counts. 
  #right now I have just counted them
  set.seed(1492)
  seq <-matrix(sample(0:1, size = 80, replace = TRUE), nc = 8)
  expect_equal(present_row(seq,seq[4,]),1)
  
  set.seed(1866)
  seq <-matrix(sample(0:1, size = 96, replace = TRUE), nc = 8)
  expect_equal(present_row(seq,seq[9,]),1)
  
  set.seed(1865)
  seq <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 8)
  expect_equal(present_row(seq,c(0,0,0,0,0,0,0,1)),0)
  
  set.seed(901)
  seq <-matrix(sample(0:1, size = 20, replace = TRUE), nc = 5)
  expect_equal(present_row(seq,c(0,0,0,0,0)),0)
  
  seq<-rbind(seq,c(1,1,1,0,1),c(1,1,1,0,1))
  expect_equal(present_row(seq,c(1,1,1,0,1)),2)
})