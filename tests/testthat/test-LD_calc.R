test_that("LD_calc works",{
  set.seed(6291)
  seq <- matrix(sample(0:1, size = 60, replace = TRUE), nc = 6)
  
  genomes = matrix2genotype(seq)
  data = genetics::LD(genomes)
  D = data$`D'`
  D_values = D[upper.tri(D)]
  
  LD_avg = mean(D_values)
  LD_max = max(D_values)
  
  r = data$`R^2`
  y = seq(2,ncol(r)-2)
  w = sapply(y, LD_w, r=r)
  w_max = max(w)
  
  df = tibble::tibble("LD_avg" = LD_avg,
                      "LD_max" = LD_max,
                      "w_max" = w_max)
  
  output = LD_calc(seq)
  expect_equal(output,df)
})