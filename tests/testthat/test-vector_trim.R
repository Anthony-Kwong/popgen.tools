#test function for vector_trim

#trim outer elements of a numeric vector

#param: x:vector
#param: cen:center of new vector
#param: k: number of elements to keep on either side of the center
#returns a trimmed vector

R_vector_trim<-function(x,cen,k){
  start = max(1,cen-k)
  end = min (length(x) , cen+k)
  y = x [start:end]
  return(y)
}

test_that("vector_trim works",{
  
  #manually checked result
  x = seq(1,10,by=1)
  #remember C++ indices start at 0
  expect_equal(x[2:6] ,vector_trim(x, cen = 4, k = 2))
  expect_equal(x[2:6], R_vector_trim(x, cen = 4, k = 2))
  
  #auto checks below
  set.seed(15)
  y = runif(0,1,n=50) %>% sort()
  cen = 23
  span = 10
  expect_equal(R_vector_trim(y, cen=cen , k=span), 
               vector_trim (y, cen=cen , k=span) )
  
  #insufficient elements on left
  set.seed(13)
  y = runif(0,1,n=100) %>% sort()
  cen = 20
  span = 30
  expect_equal(R_vector_trim(y, cen=cen , k=span), 
               vector_trim (y, cen=cen , k=span) )
  
  #insufficient elements on right
  set.seed(13)
  y = runif(0,1,n=100) %>% sort()
  cen = 85
  span = 30
  expect_equal(R_vector_trim(y, cen=cen , k=span), 
               vector_trim (y, cen=cen , k=span) )
})