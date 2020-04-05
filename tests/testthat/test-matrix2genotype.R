test_that("matrix2genotype works",{
  
  set.seed(6221941)
  SNP = 10
  nsam = 4
  seq = matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP) 
  
  g1 = genetics::genotype(seq[1,],seq[2,])
  g2 = genetics::genotype(seq[3,],seq[4,])

  output = matrix2genotype(seq)
  expect_equal(output$g_1,g1)
  expect_equal(output$g_2,g2)
  
  
  set.seed(6221242)
  SNP = 10
  nsam = 20
  seq = matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP) 
  index = seq(1, nsam, by =2)
  diploids = lapply(1:length(index), function(i) NA)
  ndip = nsam/2
  
  for ( i in 1:ndip){
    g = index[i]
    diploids[[i]] = genetics::genotype( seq[g,], seq[g+1,])
  }
  
  output = matrix2genotype(seq)
  
  for( i in 1:ndip) {
    expect_equal ( output[,i], diploids [[i]] )
  }
  
})