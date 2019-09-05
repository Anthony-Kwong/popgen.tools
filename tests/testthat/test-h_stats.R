test_that("h_stats computed correctly",{
  set.seed(1989)
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)

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
  
  count=c(1,1,1,1)+c(0,N2,N3,N4)
  freq=count/nrow(seq)
  
  output<-h_stats(seq)
  nrow(seq)
})