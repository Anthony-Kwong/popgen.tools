#Test function to compute haplotype statistics in R

#input: Binary matrix G consisting of 1's and 0's. Each column is a SNP, each row is a sample. 
#output: A vector containing h_stats (h1,h2,h12,h123)

R_hstats<-function(G){
  #compute haplotype frequencies
  hap_count<-row_freq(G)
  freq<-sort(hap_count,decreasing=TRUE)
  
  h1<-sum(freq*freq)
  h12<-h1+2*freq[1]*freq[2]
  h123<-h12+2*freq[1]*freq[3]+2*freq[2]*freq[3]
  h2<-h1-freq[1]*freq[1]
  h_stats<-c(h1,h2,h12,h123)
  
  return(h_stats)
}


test_that("hstats computed correctly",{
  set.seed(1989)
  SNP=5
  nsam=5
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:3){seq<-rbind(seq,seq[4,])}
  indices<-sample(nrow(seq))
  seq<-seq[indices, ]
  p<-c(4/8,2/8,1/8,1/8)
  h1<-sum(p*p)
  h12<-h1+2*p[1]*p[2]
  h123<-h12+2*p[1]*p[3]+2*p[2]*p[3]
  h2<-h1-p[1]^2
  ans<-c(h1,h2,h12,h123)
  expect_equal(ans,R_hstats(seq))
  expect_equal(ans,h_stats(seq))
  
  set.seed(473)
  SNP=800
  nsam=20
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:7){seq<-rbind(seq,seq[20,])}
  expect_equal(R_hstats(seq),h_stats(seq))
  
  set.seed(1688)
  SNP=1000
  nsam=30
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:7){seq<-rbind(seq[1,],seq)}
  expect_equal(R_hstats(seq),h_stats(seq))
  
  set.seed(476)
  SNP=500
  nsam=30
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:7){seq<-rbind(seq,seq[29,])}
  for(i in 1:14){seq<-rbind(seq,seq[17,])}
  indices<-sample(nrow(seq))
  seq<-seq[indices, ]
  expect_equal(R_hstats(seq),h_stats(seq))
  
  #test NA input returns NA. 
  null_win = as.numeric(NA) %>% matrix()
  check = suppressWarnings(  h_stats(null_win)  )
  null_stats = rep(NA,4) %>% as.numeric()
  expect_equal(check, null_stats)
  
  
  })

