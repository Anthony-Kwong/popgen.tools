#test function for clus_hap

R_clus_hap <- function (G, n_clus){
  uni_p = nrow(G) - sum(duplicated(G))
  
  if(uni_p < n_clus){
    n_clus = uni_p
  }
  
  #tune for the best number of clusters 
  tune_clus = factoextra::fviz_nbclust(G, kmeans, method = "silhouette",
                                       k.max = n_clus)
  tune_data = tune_clus$data
  best_clus = which.max(tune_data$y)
  
  if(best_clus < 3){
    best_clus = 3
  }
  
  G_clus = kmeans(G, centers = best_clus)
  clus_vec = G_clus$cluster
  
  return(clus_vec)
}

#add seed to kmeans

test_that("clus_hap works",{
  set.seed(1312)
  G = matrix(sample(0:1, size =5*8 , replace = TRUE), nc = 5)
  expect_equal(clus_hap(G,3), R_clus_hap(G,3))
  
  set.seed(12)
  G = matrix(sample(0:1, size =500 , replace = TRUE), nc = 10)
  expect_equal(clus_hap(G,7), R_clus_hap(G,7))
})