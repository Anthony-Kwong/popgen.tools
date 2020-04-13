test_that("LD_calc works",{
  #case 1----
  set.seed(6291)
  seq <- matrix(sample(0:1, size = 60, replace = TRUE), nc = 6)
  
  genomes = matrix2genotype(seq)
  data = genetics::LD(genomes)
  D = abs(data$`D'`)
  D_values = D[upper.tri(D)]
  
  LD_avg = mean(D_values)
  LD_max = max(D_values)
  
  r = data$`R^2`
  y = seq(2,ncol(r)-2)
  w = sapply(y, LD_w, r=r)
  w_max = max(w)
  
  r_values = r[upper.tri(r)]
  Zns = mean(r_values)
  
  df = tibble::tibble("LD_avg" = LD_avg,
                      "LD_max" = LD_max,
                      "w_max" = w_max,
                      "Zns" = Zns)
  
  output = LD_calc(seq)
  expect_equal(output,df)
  
  #check NA input returns NA
  null_win = as.numeric(NA) %>% matrix()
  check = suppressWarnings(LD_calc(null_win))
  df = tibble::tibble("LD_avg"= NA, 
                      "LD_max" = NA,
                      "w_max" = NA,
                      "Zns" = NA)
  expect_equal (check, df)
})
