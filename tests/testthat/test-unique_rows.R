test_that("fill_row function works",{
  set.seed(1492)
  seq <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 4) 
  x<-c(1,0,1,1)
  test1<-fill_row(seq,x)
  test2<-rbind(seq,x,deparse.level = 0)
  expect_equal(test1,test2)
  
  set.seed(1500)
  SNP=20
  nsam=8
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP) 
  y<-sample(c(0,1), replace=TRUE, size=SNP)
  test1<-fill_row(seq,y)
  test2<-rbind(seq,y,deparse.level = 0)
  expect_equal(test1,test2)
})

#Test function for row_count

#input: NumericMatrix A, vector x.
#return: An integer of the number of times x appeared as a row in A. 

test_count_row<-function(A,x){
  count=0
  for(i in 1:(nrow(A))){
    if(identical(A[i,],as.integer(x))==TRUE){
      count=count+1
    }
  }
  return (count)
}

test_that("row_count function works",{
  set.seed(1492)
  seq <-matrix(sample(0:1, size = 80, replace = TRUE), nc = 8)
  expect_equal(row_count(seq,seq[4,]),1)
  
  set.seed(1866)
  seq <-matrix(sample(0:1, size = 96, replace = TRUE), nc = 8)
  expect_equal(row_count(seq,seq[9,]),1)
  
  set.seed(1865)
  seq <-matrix(sample(0:1, size = 24, replace = TRUE), nc = 8)
  expect_equal(row_count(seq,c(0,0,0,0,0,0,0,1)),0)
  
  set.seed(901)
  seq <-matrix(sample(0:1, size = 20, replace = TRUE), nc = 5)
  expect_equal(row_count(seq,c(0,0,0,0,0)),0)
  
  seq<-rbind(seq,c(1,1,1,0,1),c(1,1,1,0,1))
  expect_equal(row_count(seq,c(1,1,1,0,1)),2)
  
  set.seed(1500)
  SNP=10
  nsam=50
  M <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:5){M<-rbind(M,M[1,])}
  #shuffle rows
  set.seed(1707)
  indices<-sample(nrow(M))
  M2<-M[indices, ]
  x<-M[1,]
  expect_equal(test_count_row(M2,x),row_count(M2,x))
  
})

test_that("unique_rows works",{
  set.seed(901)
  seq <-matrix(sample(0:1, size = 20, replace = TRUE), nc = 5)
  for(i in 1:6){seq<-rbind(seq,seq[4,])}
  for(j in 1:3){seq<-rbind(seq,seq[3,])}
  expect_equal(unique_rows(seq),unique(seq))
  
  set.seed(1836)
  SNP=100
  nsam=20
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:10){seq<-rbind(seq,seq[2,])}
  for(j in 1:5){seq<-rbind(seq,seq[18,])}
  expect_equal(unique_rows(seq),unique(seq))
  
  # set.seed(800)
  # seq <-matrix(sample(0:1, size = 16, replace = TRUE), nc = 4)
  # seq<-rbind(seq,seq[2,])
  # seq<-rbind(seq,seq[4,])
  # seq<-rbind(seq,seq[3,])
  # expect_equal(unique_rows(seq),c(1,2,2,2))
})

#Test function for row_freq

#Finds the frequencies of each unique row in a given matrix

#input: NumericMatrix A
#return: a vector of the frequencies of each unique row

test_row_freq<-function(A){
  B<-unique(A)
  num_hap=nrow(B)
  freq=rep(NA,num_hap)
  
  for(i in 1:num_hap){
    freq[i]=row_count(A,B[i,])
  }
  
  freq=freq/nrow(A)
  
  return(freq)
}

test_that("row_freq function works",{
  set.seed(1836)
  SNP=5
  nsam=5
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:3){seq<-rbind(seq,seq[1,])}
  expect_equal(row_freq(seq),test_row_freq(seq))
  
  set.seed(1066)
  SNP=500
  nsam=50
  seq <-matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP)
  for(i in 1:10){seq<-rbind(seq,seq[7,])}
  expect_equal(row_freq(seq),test_row_freq(seq))
})

test_that("vec_sort works",{
  set.seed(1187)
  x<- sample(1:35, 10, replace=FALSE)
  expect_equal(sort(x,decreasing = T),vec_sort(x))
  
  y<-runif(n=50,min=0,max=1)
  expect_equal(sort(y,decreasing = T),vec_sort(y))
})

