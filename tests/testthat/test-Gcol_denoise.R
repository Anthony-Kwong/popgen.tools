#test function for Gcol_denoise

#this version uses 1 indexing unlike the Cpp version

R_Gcol_denoise <- function(G,col){
  pwd = Gcol_flip(G,col-1)
  flip_index = which.min(pwd)
  el = G[flip_index,col]
  if(el == 0){
    G[flip_index,col] = 1
  } else if(el == 1){
    G[flip_index,col] = 0
  } else {
    stop("Invalid genome matrix")
  }
  return(G)
}

test_that("Gcol_denoise works",{
  set.seed(1232)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  ans = Gcol_denoise(G,0)
  R_says = R_Gcol_denoise(G,1)
  G[3,1] = 1
  expect_equal(G, ans)
  expect_equal(ans, R_says)
  
  set.seed(122)
  G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
  ans = Gcol_denoise(G,0)
  R_says = R_Gcol_denoise(G,1)
  G[3,1] = 1
  expect_equal(G, ans)
  expect_equal(R_says, ans)
  
  set.seed(121232)
  G = matrix(sample(0:1, size =180 , replace = TRUE), nc = 12)
  expect_equal(Gcol_denoise(G,11),R_Gcol_denoise(G,12))
})
