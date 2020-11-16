#test function for clus_hstats
R_clush <- function(x){
  num_hap = max(x)
  hap_count = sapply(seq(1:num_hap), function(i){count(x,i)})
  hap_freq = hap_count/length(x)
  hap_freq = sort(hap_freq, decreasing = T)
  
  #need to sort frequencies
  hap_sq = hap_freq*hap_freq
  h1 = sum(hap_sq)
  h2 = sum(hap_sq[2:num_hap])
  
  #account for having fewer clusters to ensure R doesn't read the hap_sq vector backwards
  if(length(hap_sq) < 3){
    rem = 0
  } else {
    rem = hap_sq[3:num_hap]
  }
  h12 = (hap_freq[1]+hap_freq[2])^2 + sum(rem)
  
  if(length(hap_sq) < 4){
    rem = 0
  } else {
    rem = hap_sq[4:num_hap]
  }
  
  h123 = (hap_freq[1]+hap_freq[2]+hap_freq[3])^2 + sum(rem)
  
  ans = c(h1,h2,h12,h123)
  return(ans)
}

test_that("clus_hstats works",{
  #manual check
  set.seed(131211)
  G = matrix(sample(0:1, size =10*8 , replace = TRUE), nc = 10)
  x = clus_hap(G,4)
  h1 = 2*(3/8)^2 + (2/8)^2
  h2 = h1 - (3/8)^2
  h12 = (6/8)^2 + (2/8)^2
  h123 = 1
  expect_equal(clus_hstats(x), c(h1, h2, h12, h123))
  expect_equal(clus_hstats(x), R_clush(x))
  
  set.seed(1311)
  G = matrix(sample(0:1, size =10 * 10 , replace = TRUE), nc = 10)
  x = clus_hap(G,3)
  h1 = 2*(4/10)^2 + (2/10)^2
  h2 = h1 - (4/10)^2
  h12 = (8/10)^2 + (2/10)^2
  h123 = 1
  expect_equal(clus_hstats(x), c(h1, h2, h12, h123))
  expect_equal(clus_hstats(x), R_clush(x))
  
  set.seed(11)
  G = matrix(sample(0:1, size =10 * 10 , replace = TRUE), nc = 10)
  x = clus_hap(G,6)
  h1 = (4/10)^2 + 3*(2/10)^2
  h2 = h1 - (4/10)^2
  h12 = (6/10)^2 + 2*(2/10)^2
  h123 = (8/10)^2 + (2/10)^2
  expect_equal(clus_hstats(x), c(h1, h2, h12, h123))
  expect_equal(clus_hstats(x), R_clush(x))
  
  #autocheck
  set.seed(131211)
  G = matrix(sample(0:1, size =50*20 , replace = TRUE), nc = 50)
  x = clus_hap(G,10)
  expect_equal(clus_hstats(x), R_clush(x))
  
  
  set.seed(52)
  G = matrix(sample(0:1, size =30*25 , replace = TRUE), nc = 30)
  x = clus_hap(G,15)
  expect_equal(clus_hstats(x), R_clush(x))
})