#R test function for finding the index for the vector element closest to target value

#input: x vector to search through
#output: target some target value
#return: index of the closest vector element -1. 
#(to account for indices starting at 0 in CPP)

R_find_index<-function(x,target){
  v<-abs(x-target)
  index<-which.min(v)-1
  return(index)
}

test_that("find_index works",{
  #manual check
  set.seed(1444)
  x<-runif(0,1,n=5) %>% sort()
  d<-0.2
  answer<-0
  expect_equal(R_find_index(x,d),answer)
  expect_equal(find_index(x,d),answer)
  
  set.seed(144)
  x<-runif(0,1,n=10) %>% sort()
  d<-0.5
  answer<-5-1
  expect_equal(R_find_index(x,d),answer)
  expect_equal(find_index(x,d),answer)
  
  #auto check
  set.seed(52)
  x<-runif(0,1,n=550) %>% sort()
  d<-0.8
  expect_equal(R_find_index(x,d),find_index(x,d))
})