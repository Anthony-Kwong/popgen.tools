test_that("norm_vec works",{
  e<-5
  std=1
  set.seed(1066)
  x<-rnorm(100,mean=e,sd=std)
  output<-norm_vec(x)
  actual<-(x-mean(x))/(var(x)^0.5)
  expect_equal(output,actual)
})