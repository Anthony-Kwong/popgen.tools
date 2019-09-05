test_that ("simulation object constructly correctly", {
  cmd="~/command"
  seeds=c(1,2)
  segsites=10
  positions=c(0.4,0.6)
  set.seed(2019)
  genome_matrix=matrix(sample(0:1, size = 9, replace = TRUE), nc = 3)
  sweep="hard"
  select_coeff=0.1
  sim=sim_obj(cmd,seeds,segsites,positions,genome_matrix,sweep,select_coeff)

  expect_equal(cmd,sim$cmd)
  expect_equal(seeds,sim$seeds)
  expect_equal(segsites,sim$num_seg)
  expect_equal(positions,sim$pos)
  expect_equal(genome_matrix,sim$genomes)
  expect_equal(sweep,sim$sweep)
  expect_equal(select_coeff,sim$s)
})


