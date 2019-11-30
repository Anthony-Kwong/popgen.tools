test_that("string_labels works",{
  out<-string_labels("foo",3)
  actual<-c("foo_1","foo_2","foo_3")
  expect_equal(out,actual)
  
  str<-"smoo"
  num<-30
  
  out<-string_labels(str,num)
  actual<-stringr::str_c(str,1:num,sep="_")
  expect_equal(out,actual)
})