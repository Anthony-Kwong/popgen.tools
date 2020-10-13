test_that("Gcol_denoise works",{
  set.seed(1232)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  Gcol_flip(G,0)
  Gcol_denoise(G,0)
})