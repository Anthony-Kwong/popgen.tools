#test function for Gcol_flip

#this version uses 1 indexing unlike the Cpp version

R_Gcol_flip <- function(G,col){
  pd = rep(NA,nrow(G))
  for(i in 1:nrow(G)){
    
    K = G
    el = K[i,col]
   
    if(el == 0){
      K[i,col] = 1
    } else if (el == 1){
      K[i,col] = 0
    } else {
      stop("Input matrix has an element that is neither 0 nor 1.")
    }
    
    pd[i] = pairwise_diff(K)
  }
  return(pd)
}

test_that("Gcol_flip works",{
  set.seed(12)
  G = matrix(sample(0:1, size =3*3 , replace = TRUE), nc = 3)
  expect_equal(Gcol_flip(G,0),c(6,6,6))
  expect_equal(R_Gcol_flip(G,1), Gcol_flip(G,0))
  
  
  set.seed(12)
  G = matrix(sample(0:1, size =20 , replace = TRUE), nc = 5)
  expect_equal(Gcol_flip(G,2), c(14,14,14,10))
  expect_equal(R_Gcol_flip(G,3), Gcol_flip(G,2))
  
  set.seed(112)
  G = matrix(sample(0:1, size =80 , replace = TRUE), nc = 8)
  expect_equal(Gcol_flip(G,7), R_Gcol_flip(G,8))
  
  
  set.seed(11223)
  G = matrix(sample(0:1, size =120 , replace = TRUE), nc = 12)
  expect_equal(Gcol_flip(G,0), R_Gcol_flip(G,1))
})