test_that("Gcol_flip works",{
  set.seed(12)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  expect_equal(Gcol_flip(G,0),c(6,6,6))
  
  set.seed(12)
  G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
  expect_equal(Gcol_flip(G,2), c(14,14,14,10))
})