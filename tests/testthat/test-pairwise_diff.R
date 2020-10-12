test_that("pairwise_diff works",{
  set.seed(42)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  expect_equal(pairwise_diff(G),4)
  
  set.seed(54)
  G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
  expect_equal(pairwise_diff(G),17)
})