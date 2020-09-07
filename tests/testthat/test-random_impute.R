#defunction R test function because I couldn't get the random seeds to align with Rcpp

# R_random_impute <- function (G){
#   output = G
#   sites = ncol(G)
#   for(c in 1:sites){
#     col = G[,c]
#     na_pos = which(is.na(col))
#     num_na = length(na_pos)
#     col = purrr::discard(col,is.na)
#     for (r in 1:num_na){
#       imp = runif(1,min=0,max=length(col))
#       if(num_na==0){
#         break
#       }
#       print(imp)
#       imp_index = round(imp) #account for 1 indexing for R
#       #print(imp_index)
#       output[na_pos[r],c] = col[imp_index + 1]
#     }
#   }
#   return (output)
# }

test_that("random_impute works",{
  #manual check, got the random imputations using Rcout statements
  set.seed(2)
  G = matrix(sample(0:1, size = 25, replace = TRUE), nc = 5)
  G[1,1] = NA
  G[3,4] = NA
  G[5,2] = NA
  set.seed(5)
  out = random_impute(G)
  
  G[1,1] = 0
  G[5,2] = 0
  G[3,4] = 1
  expect_equal(out, G)
  
  #manual check 2
  set.seed(22)
  G = matrix(sample(0:1, size = 20, replace = TRUE), nc = 5)
  G = age_DNA(G, missing_rate = 0.25, seed = 20)
  
  set.seed(3)
  out = random_impute(G)
  G[,3] = rep(0,4)

  G[3,4] = 0
  G[4,4] = 1
  
  G[2,5] = 0
  G[3,5] = 1
  expect_equal(G,out)
  
  #manual check with a column of NAs, make sure it is convereted to 0's. 
  set.seed(22)
  G = matrix(sample(0:1, size = 20, replace = TRUE), nc = 5)
  G[,4] = rep(NA,4)
  
  output = random_impute(G)
  G[,4] = rep(0,4)
  expect_equal(output,G)
  
})