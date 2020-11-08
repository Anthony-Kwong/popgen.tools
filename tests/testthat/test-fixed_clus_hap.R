#test function for fixed_clus_hap

R_clus <- function(G, k, seed = NA) {
  
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1)
  }
  
  uni_dp <- unique_rows(G) %>% nrow()
  
  if(uni_dp < k){
    k = uni_dp
  }
  
  set.seed(seed)
  G_clus <- kmeans(G, centers = k)
  clus_vec <- G_clus$cluster
  
  return(clus_vec)
}

test_that("fixed_clus_hap function works",{
  set.seed(1311)
  G = matrix(sample(0:1, size =5*8 , replace = TRUE), nc = 5)
  expect_equal(fixed_clus_hap(G, 3, seed = 32), R_clus(G, 3, seed = 32))
  
  set.seed(123)
  G = matrix(sample(0:1, size =500 , replace = TRUE), nc = 10)
  expect_equal(fixed_clus_hap(G,8, seed = 12), R_clus(G,8, seed = 12))
})