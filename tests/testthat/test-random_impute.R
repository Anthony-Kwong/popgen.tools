R_random_impute <- function (G){
  output = G
  sites = ncol(G)
  nsam = nrow(G)
  for(c in 1:sites){
    col = G[,c]
    na_pos = which(is.na(col))
    num_na = length(na_pos)
    col = purrr::discard(col,is.na)
    for (r in 1:num_na){
      imp = runif(1,min=0,max=length(col))
      imp_index = round(imp)
      output[na_pos[r],c] = col[imp_index]
    }
  }
  return (output)
}

test_that("random_impute works",{
  #manual check
  set.seed(2)
  G = matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
  G[1,1] = NA
  G[3,4] = NA
  G[5,2] = NA
  set.seed(3)
  out = random_impute(G)
  set.seed(3)
  R_says = R_random_impute(G)
  
  G[1,1] = 0
  G[5,2] = 0
  G[3,4] = 1
  expect_equal(out, G)
  expect_equal(out, R_says)
})