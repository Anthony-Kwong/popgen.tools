test_that("sub_win works",{
  
  set.seed(7)
  SNP=11
  num_windows=3
  seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
  output=sub_win(G=seq,num_windows=num_windows)
  width=ceiling(SNP/num_windows)
  
  truth=list()
  for(i in 1:num_windows){
    if(i<num_windows){
      truth[[i]]=seq[,((i-1)*width+1):(i*width)]
    } else {
      truth[[i]]=seq[,((i-1)*width+1):ncol(seq)]
    }
  }
  
  #checking the 2 lists are equal
  for(i in 1:num_windows){
    expect_equal(all.equal(output[[i]],truth[[i]]),T)
   }
})
  