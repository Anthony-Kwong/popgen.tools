test_that("Converted away from scientific notation correctly",{
  x<-2.4e8
  test<-no_scientific(x)
  expect_equal(test,format(x,scientific = F,trim = T))

  x<-0
  test<-no_scientific(x)
  expect_equal(test,format(x,scientific = F,trim = T))

  x<--6.7e-10
  test<-no_scientific(x)
  expect_equal(test,format(x,scientific = F,trim = T))
} )

