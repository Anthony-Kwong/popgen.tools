#test function for clus_hstats
R_clush <- function(x){
  num_hap = max(x)
  hap_count = sapply(seq(1:num_hap), function(i){count(x,i)})
  hap_freq = hap_count/length(x)
  
  hap_sq = hap_freq*hap_freq
  h1 = sum(hap_sq)
  h2 = sum(hap_sq[2:num_hap])
  
  #remove NA values for cases where there is <3 clusters
  if(length(hap_sq) < 3){
    rem = 0
  } else {
    rem = hap_sq[3:num_hap]
    #rem = rem[!is.na(rem)]
  }
  h12 = (hap_freq[1]+hap_freq[2])^2 + sum(rem)
  
  if(length(hap_sq) < 3){
    rem = 0
  } else {
    rem = hap_sq[4:num_hap]
  }
  h123 = (hap_freq[1]+hap_freq[2]+hap_freq[3])^2 + sum(rem)
  
  ans = c(h1,h2,h12,h123)
  return(ans)
}

test_that("clus_hstats works",{
  set.seed(131211)
  G = matrix(sample(0:1, size =10*8 , replace = TRUE), nc = 10)
  x = clus_hap(G,6)
  expect_equal(clus_hstats(x), c(0.5, 0.25, 1, NA))
  expect_equal(clus_hstats(x), R_clush(x))
  
  set.seed(131211)
  G = matrix(sample(0:1, size =50*20 , replace = TRUE), nc = 50)
  x = clus_hap(G,10)
  expect_equal(clus_hstats(x), R_clush(x))
})