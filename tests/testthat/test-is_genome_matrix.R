test_that("is_genome_matrix works",{
  x = "2"
  expect_equal(is_genome_matrix(x), FALSE)
  
  x = rbind(c(T,T), c(F,F)) 
  expect_equal(is_genome_matrix(x), FALSE)
  
  x = rbind (c(1,0,0,0), c(0,1,0,2), c(0,0,1,0))
  expect_equal(is_genome_matrix(x), FALSE)
  
  x = rbind (c(1,0,0,0), c(0,1,0,1), c(0,0,0,0))
  expect_equal(is_genome_matrix(x), FALSE)
  
  x = rbind (c(1,0,0,1,0), c(0,1,1,1,0), c(0,0,0,0,1))
  expect_equal(is_genome_matrix(x), TRUE)
  
  set.seed(1444)
  x <-matrix(sample(0:1, size =20 , replace = TRUE), nc = 4)
  expect_equal(is_genome_matrix(x), TRUE)
  
  x <- cbind (x, rep(1,5))
  expect_equal(is_genome_matrix(x), F)
  
})