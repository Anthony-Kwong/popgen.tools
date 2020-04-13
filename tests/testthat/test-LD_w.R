test_that("LD_w works",{
  
  #manual check
  set.seed(2)
  mat <- matrix(runif(16), 4, 4)
  mat[upper.tri(mat)==F] <- NA
  nsam = ncol(mat)
  i = 2
  L = mat[,1:i]
  R = mat[, (i+1):nsam]
  
  top_s = sum(L,na.rm = T) + sum(R[(i+1):nsam,],na.rm = T)
  top_n = choose(i,2) + choose (nsam - i, 2)
  top = top_s / top_n
  
  bot = sum(R[1:i,]) / ( i*(nsam-i) )
  w = top/bot
  
  expect_equal(w, LD_w(mat,i))
})
