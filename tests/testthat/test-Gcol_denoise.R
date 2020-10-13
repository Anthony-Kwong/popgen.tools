test_that("Gcol_denoise works",{
  set.seed(1232)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  ans = Gcol_denoise(G,0)
  G[3,1] = 1
  expect_equal(G, ans)
  
  set.seed(122)
  G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
  ans = Gcol_denoise(G,0)
  G[1,3] = 1
  expect_equal(G, Gcol_denoise(G,2))
})