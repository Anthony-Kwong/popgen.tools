#test function for G_flip

R_G_flip <- function(G){
  snp = ncol(G)
  for(i in 1:snp){
    G = Gcol_denoise(G,i-1)
  }
  return(G)
}

test_that("G_flip function works",{
  set.seed(1232)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  expect_equal(G_flip(G),cbind(c(1,1,1),c(0,0,0),c(1,1,1)))
  expect_equal(G_flip(G), R_G_flip(G))
  
  set.seed(132)
  G = matrix(sample(0:1, size =15 , replace = TRUE), nc = 5)
  expect_equal(G_flip(G),rbind( c(0,1,0,0,1),c(0,1,0,0,1), c(0,1,0,0,1)) )
  expect_equal(G_flip(G), R_G_flip(G))
  
  set.seed(1322)
  G = matrix(sample(0:1, size =40*15 , replace = TRUE), nc = 15)
  expect_equal(G_flip(G), R_G_flip(G))
  
})