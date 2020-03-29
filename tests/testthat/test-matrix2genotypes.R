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
})