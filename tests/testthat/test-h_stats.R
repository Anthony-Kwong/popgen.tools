test_that("h_stats computed correctly",{
  set.seed(1989)
  #default, this generates 4 unique rows
  seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)
  #might turn this bit into a function later
  copy2<-6
  copy3<-4
  copy4<-5
  
  for(i in 1:copy2){
    seq<-rbind(seq,seq[2,])
  }
  
  seq<-rbind(seq,seq[2,])
  seq<-rbind(seq,seq[4,])
  seq<-rbind(seq,seq[3,])
})