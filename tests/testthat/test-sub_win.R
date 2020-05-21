#test function for win_split (R version).
#
# When ncol(seq) is not divisible by num_windows, the remaining columns are sent to the last subwindow.
#
#param: seq: a numeric matrix, num_windows: the number of sub windows to split matrix
#return: A list of 2 elements. 1 is a list of subwindows that make up the original matrix.
# 2 is a vector of column indices, indicating the boundary points for each subwindow. 

R_win_split<-function(seq,num_windows){
  wins=list()
  width=(ncol(seq)/num_windows) %>% floor()
  col_index = rep(NA, num_windows + 1 )
  col_index[1] = 1
  for(i in 1:num_windows){
    if(i<num_windows){
      wins[[i]]=seq[,((i-1)*width+1):(i*width)]
      col_index[i+1] = i*width
    } else {
      wins[[i]]=seq[,((i-1)*width+1):ncol(seq)]
      col_index[i+1] = ncol(seq)
    }
  }
  output = list (windows = wins,
                 bounds = col_index)
  return (output)
}

test_that("sub_win works",{
  
  #manual check
  set.seed(7)
  SNP=11
  num_windows=3
  seq <-matrix(sample(0:1, size = SNP*3, replace = TRUE), nc = SNP)
  output=sub_win( G = seq, num_windows = num_windows)
  R_says = R_win_split( seq, num_windows)
  
  expect_equal(R_says, output)
  
  #for cases where the number of subwindows requested is bigger than the number of SNPs. 
  set.seed(2261941 )
  SNP=8
  num_windows=10
  seq = matrix(sample(0:1, size = SNP*20, replace = TRUE), nc = SNP)
  output = sub_win(G = seq, num_windows = num_windows)
  expect_equal(output, list())
  
  set.seed(2261941)
  SNP=80
  num_windows=10
  seq <-matrix(sample(0:1, size = SNP*4, replace = TRUE), nc = SNP)
  output=sub_win(G=seq,num_windows=num_windows)
  R_says<-R_win_split(seq,num_windows)
  expect_equal(output, R_says)
  
})
  