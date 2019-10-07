#test function for win_split (R version).

#param: seq: a numeric matrix, num_windows: the number of sub windows to split matrix
#return: a list of subwindows
R_win_split<-function(seq,num_windows){
  wins=list()
  width=(ncol(seq)/num_windows) %>% ceiling()
  for(i in 1:num_windows){
    if(i<num_windows){
      wins[[i]]=seq[,((i-1)*width+1):(i*width)]
    } else {
      wins[[i]]=seq[,((i-1)*width+1):ncol(seq)]
    }
  }
  return (wins)
}

test_that("sub_win works",{
  
  set.seed(7)
  SNP=11
  num_windows=3
  seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
  output=sub_win(G=seq,num_windows=num_windows)
  #width=ceiling(SNP/num_windows)
  
  truth<-R_win_split(seq,num_windows)
  
  #checking the 2 lists are equal
  for(i in 1:num_windows){
    expect_equal(all.equal(output[[i]],truth[[i]]),T)
  }
  
  #for cases where the number of subwindows requested is bigger than the number of SNPs. 
  set.seed(2261941)
  SNP=8
  num_windows=10
  seq <-matrix(sample(0:1, size = SNP*20, replace = TRUE), nc = SNP)
  output=sub_win(G=seq,num_windows=num_windows)
  expect_equal(output,list())
})
  