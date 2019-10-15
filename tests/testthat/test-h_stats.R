#create(h1,h2,h12,h123)
R_hstats<-function(G){
  H<-unique(G)
  freq<-rep(NA,nrow(H))
  
  freq<-freq*nrow(G)
  
  #compute h stats
  h1<-sum(freq*freq)
  h12<-h12+2*freq[1]*freq[2]
  h123<-h
}
#c1 <- which(M[, 1] == v[1])


test_that("h_stats computed correctly",{
  set.seed(1989)
  SNP=4
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = SNP)

  #we will duplicate the following rows
  copy2<-2
  copy3<-3
  copy4<-4
  
  N2<-3
  N3<-5
  N4<-7
  
  seq<-row_replicate(seq,row=copy2,n=N2)
  seq<-row_replicate(seq,row=copy3,n=N3)
  seq<-row_replicate(seq,row=copy4,n=N4)
  
  count=rep(1,SNP)
  count[copy2]=count[copy2]+N2
  count[copy3]=count[copy3]+N3
  count[copy4]=count[copy4]+N4
  
  freq=count/nrow(seq)
  freq=sort(freq,decreasing = TRUE)
  
  h1=freq %*% freq
  h2=h1-freq[1]^2
  h12=h1+2*freq[1]*freq[2]
  h123=h12+2*freq[1]*freq[3]+2*freq[2]*freq[3]
  
  
  h=c(h1,h2,h12,h123)
  output<-h_stats(seq)
  
  expect_equal(output,h)
})

test_that("h_stats computed correctly",{
  set.seed(1929)
  #default, this generates 5 unique rows. This is important for testing. Code breaks otherwise. 
  SNP=5
  seq <-matrix(sample(0:1, size = 25 , replace = TRUE), nc = SNP)
  
  #we will duplicate the following rows
  copy2<-3
  copy3<-2
  copy4<-1
  
  N2<-8
  N3<-10
  N4<-7
  
  seq<-row_replicate(seq,row=copy2,n=N2)
  seq<-row_replicate(seq,row=copy3,n=N3)
  seq<-row_replicate(seq,row=copy4,n=N4)
  
  count=rep(1,SNP)
  count[copy2]=count[copy2]+N2
  count[copy3]=count[copy3]+N3
  count[copy4]=count[copy4]+N4
  
  freq=count/nrow(seq)
  freq=sort(freq,decreasing = TRUE)
  
  h1=freq %*% freq
  h2=h1-freq[1]^2
  h12=h1+2*freq[1]*freq[2]
  h123=h12+2*freq[1]*freq[3]+2*freq[2]*freq[3]
  
  
  h=c(h1,h2,h12,h123)
  output<-h_stats(seq)
  
  expect_equal(output,h)
})