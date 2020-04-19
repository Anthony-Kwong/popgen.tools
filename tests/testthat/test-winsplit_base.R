#test function for winsplit_base
#
#Breaks a genome matrix into equal chunks based on genome length
#
#param: NumericMatrix G, binary genome matrix
#param: NumericVector pos, sorted vector with values between 0,1
#param: integer n. number of chunks to make.
#return: A list of matrices.

R_winsplit_base<-function(G, pos, n){
  chrom_len = tail(pos,1)-pos[1]
  len = chrom_len/n
  start_pos=pos[1]
  start_indices=rep(NA,n+1)
  start_indices[1]=1
  
  for(i in 2:(n+1)){
    target = start_pos + (i-1)*len
    dist = pos - target
    index = which(dist >= 0)[1] #find the first index beyond boundary
    start_indices[i] = index
  }
  
  wins = list()
  start= start_indices[1]

  for(j in 1:n){
    if(j == n){
      end = start_indices[j+1]
    } else {
      end=start_indices[j+1]-1
    }
    
    #return NA if a block has no SNPs
    if(start_indices[j] == start_indices[j+1]){
      wins[[j]]= as.numeric(NA) %>% matrix()
      next
    }
    
    wins[[j]]=G[,start:end] %>% as.matrix()
    start = end + 1
  }
  
  return(wins)
}


test_that("winsplit_base works",{
  
  #manual check, no trim ----
  set.seed(20)
  SNP=10
  seq <-matrix(sample(0:1, size = SNP*6, replace = TRUE), nc = SNP)
  pos <- runif(0,1,n=SNP) %>% sort()
  pos[1] = 0
  pos[10] = 1
  C_says <- winsplit_base(seq,pos,n=3)
  R_says <- R_winsplit_base(seq,pos,n=3)
  
  ans=list()
  ans[[1]] = seq[,1:4]
  ans[[2]] = seq[,5:7]
  ans[[3]] = seq[,8:10]
  
  expect_equal(C_says, ans)
  expect_equal(R_says, R_says)
  
  #test for trimmed simulations, manual ----
  set.seed(21)
  SNP=20
  nwins=3
  nsam=5
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  pos <- runif(0,1,n=SNP) %>% sort()
  cen = 10
  k = 5
  
  G = window_trim(seq,cen,k)
  pos_trim = vector_trim(pos, cen, k)
  C_says = winsplit_base(G, pos_trim, nwins)
  R_says = R_winsplit_base( G, pos_trim, nwins)

  ans = list()
  ans[[1]] = G[,1:4]
  ans[[2]] = G[,5] %>% as.matrix()
  ans[[3]] = G[,6:11]
  
  expect_equal(R_says, ans)
  
  #expect_equal(C_says, R_says)
  
  #manual check, no SNPs in a block----
  set.seed(21)
  SNP = 20
  nwins = 4
  nsam = 5
  seq <- matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  pos1 <- runif(0,0.25,n=(SNP/2)) %>% sort()
  pos2 <- runif(0.75,1,n=(SNP/2)) %>% sort()
  pos = c(pos1, pos2) 
  
  ans = list()
  ans[[1]] = seq[,1:10]
  ans[[2]] = NA %>% as.numeric %>% as.matrix()
  ans[[3]] = NA %>% as.numeric %>% as.matrix()
  ans[[4]] = seq[,11:20]
  
  #C_says = suppressWarnings ( winsplit_base( seq, pos, nwins) )
  R_says = R_winsplit_base( seq, pos, nwins)
  
  expect_equal(ans,R_says)
  
  #autocheck, untrimmed
  
  set.seed(21)
  SNP=250
  nwins=10
  seq <-matrix(sample(0:1, size = SNP*50, replace = TRUE), nc = SNP)
  pos <- runif(0,1,n=SNP) %>% sort()
  C_says <- winsplit_base(seq,pos,n=nwins)
  R_says <- R_winsplit_base(seq,pos,n=nwins)
  
  expect_equal(C_says, R_says)
  
  

  
  #auto check for trimmed sims
  set.seed(21)
  SNP = 500
  nwins = 10
  nsam = 30
  seq <- matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  pos <- runif(0,1,n=SNP) %>% sort()
  cen = 250
  k = 125
  
  G = window_trim(seq,cen,k)
  pos_trim = vector_trim(pos, cen, k)
  C_says = winsplit_base(G, pos_trim, nwins)
  R_says = R_winsplit_base( G, pos_trim, nwins)
  
  expect_equal(C_says, R_says)
  
  #checking effect of having no SNPs in a block
  
  #auto check for trimmed sims
  set.seed(21)
  SNP = 20
  nwins = 4
  nsam = 5
  seq <- matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  pos1 <- runif(0,0.25,n=(SNP/2)) %>% sort()
  pos2 <- runif(0.75,1,n=(SNP/2)) %>% sort()
  pos = c(pos1, pos2) 
  
  C_says = suppressWarnings ( winsplit_base( seq, pos, nwins) )
  R_says = R_winsplit_base( seq, pos, nwins)
  
  expect_equal(C_says, R_says)
  
})
