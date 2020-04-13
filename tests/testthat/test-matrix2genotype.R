test_that("matrix2genotype works",{
  
  #case 1----
  set.seed(6221941)
  SNP = 10
  nsam = 6
  seq = matrix(sample(0:1, size = SNP*nsam, replace = TRUE), nc = SNP) 
  
  g1 = genetics::genotype(seq[1,],seq[2,])
  g2 = genetics::genotype(seq[3,],seq[4,])
  g3 = genetics::genotype(seq[5,],seq[6,])
  

  output = matrix2genotype(seq)
  df = data.frame(g1,g2,g3) %>% t() %>% genetics::makeGenotypes()
  expect_equal(df,output)

  #case 2 ----
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
  gen_df = do.call(data.frame, diploids)
  names(gen_df) = paste0("g",1:length(diploids))
  gen_df = genetics::makeGenotypes(t(gen_df))
  
  expect_equal(gen_df,output)
  
  #check NA input returns NA
  check = suppressWarnings(matrix2genotype(NA))
  expect_equal (check, NA)
})