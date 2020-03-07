test_that ("simulation object constructly correctly", {
  
  #constant pop size
  cmd="~/command"
  seeds=c(1,2)
  segsites=10
  positions=c(0.4,0.6)
  set.seed(2019)
  genome_matrix=matrix(sample(0:1, size = 9, replace = TRUE), nc = 3)
  sweep="hard"
  select_coeff=0.1
  fix_time=3
  sim=sim_obj(cmd,seeds,segsites,positions,genome_matrix,sweep,select_coeff,fix_time=fix_time)

  expect_equal(cmd,sim$cmd)
  expect_equal(seeds,sim$seeds)
  expect_equal(segsites,sim$num_seg)
  expect_equal(positions,sim$pos)
  expect_equal(genome_matrix,sim$genomes)
  expect_equal(sweep,sim$sweep)
  expect_equal(select_coeff,sim$s)
  expect_equal(fix_time,sim$fix_time)
  
  #check default bottleneck stats are recorded
  expect_equal(0,sim$bottle_time1)
  expect_equal(0,sim$bottle_time2)
  expect_equal(1,sim$bottle_size1)
  expect_equal(1,sim$bottle_size2)
  
})

test_that("simuation object constructed correctly", {
  
  #bottleneck scenario
  cmd="~/command"
  seeds=c(1,2)
  segsites=10
  positions=runif(0,1,n=segsites) %>% sort()
  set.seed(2019)
  genome_matrix=matrix(sample(0:1, size = 90, replace = TRUE), nc = 10)
  sweep="neutral"
  select_coeff=0
  fix_time=5
  bt1=5
  bs1=0.4
  bt2=10
  bs2=0.8
  
  sim=sim_obj(cmd,seeds,segsites,positions,genome_matrix,sweep,
              select_coeff,fix_time=fix_time,bottle_time1 = bt1, 
              bottle_size1 = bs1,bottle_time2 = bt2, bottle_size2 = bs2)
  
  expect_equal(cmd,sim$cmd)
  expect_equal(seeds,sim$seeds)
  expect_equal(segsites,sim$num_seg)
  expect_equal(positions,sim$pos)
  expect_equal(genome_matrix,sim$genomes)
  expect_equal(sweep,sim$sweep)
  expect_equal(select_coeff,sim$s)
  expect_equal(fix_time,sim$fix_time)
  
  #check default bottleneck stats are recorded
  expect_equal(bt1,sim$bottle_time1)
  expect_equal(bt2,sim$bottle_time2)
  expect_equal(bs1,sim$bottle_size1)
  expect_equal(bs2,sim$bottle_size2)
})


