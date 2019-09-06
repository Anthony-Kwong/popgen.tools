library(stringr)

test_that("Command entered correctly", {

  #input parameters
  mu=1e-8
  recomb_rate=1e-9
  Ne=1000000
  nSites=1e5
  samplesize=20
  s=0.1
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
  cmd<-sim$cmd

  alpha=no_scientific(2*Ne*s) #scaled strength of selection
  theta=no_scientific(4*Ne*mu*nSites) #scaled mutation rate
  rho=no_scientific(4*Ne*recomb_rate)  #recomb_rate is the probability of a cross over per basepair of sequence being modelled.
  tau= no_scientific(fix/(4*Ne)) #scaled time for fixation

  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(nSites),"-t",theta,"-r", rho,
                 "-a", alpha, "-ws", tau)
  expect_equal(cmd,test_cmd)

  seeds=c(1,2)
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  cmd<-sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(nSites),"-t",theta,"-r", rho,"-d", seeds[1],seeds[2],
                 "-a", alpha, "-ws", tau)
  expect_equal(cmd,test_cmd)

  #testing the neutral sweep command

  sweep="neutral"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
  input_cmd=sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(nSites),"-t",theta,"-r", rho,
                 "-wn", tau)
  expect_equal(input_cmd,test_cmd)
})

test_that("Seed extraction successful",{

  #input parameters
  mu=runif(1,1e-8,1e-7)
  recomb_rate=runif(1,1e-10,1e-9)
  Ne=runif(1,1000,1000000)
  nSites=1e5
  samplesize=runif(1,1,50)
  s=runif(1,1e-3,0.1)
  fix=runif(1,1,5)
  discoal_path="~/work/programs/discoal/discoal"
  seed1=floor(runif(1,1,9e6))
  seed2=floor(runif(1,1,9e6))
  seeds=c(seed1,seed2)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  used_seed<-sim$seeds
  expect_equal(used_seed,seeds)
})

test_that("Segsites extraction successful",{
  #input parameters
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000000
  nSites=1e5
  samplesize=15
  s=0.1
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  seeds=c(2019,1688)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  extracted_seg=sim$num_seg
  print(extracted_seg)
  #Ran in terminal to get this
  actual_seg=9
  expect_equal(extracted_seg,actual_seg)

  samplesize=1
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  expect_equal(sim,"Simulation produced no segregating sites")

})

test_that("Postions extraction successful",{
  #input parameters
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000000
  nSites=1e5
  samplesize=15
  s=0.1
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  seeds=c(1453,1688)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  extracted_pos=sim$pos
  #Ran the command in terminal to get these numbers
  actual_pos=c(0.074125, 0.462233, 0.539124, 0.572956, 0.608014, 0.661452, 0.673324, 0.715790, 0.896508 )
  expect_equal(extracted_pos,actual_pos)
})

test_that("Genome matrix extraction successful",{
  #input parameters
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000000
  nSites=1e5
  samplesize=5
  s=0.1
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  seeds=c(1688,1707)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,nSites=nSites,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)

  #Obtained via directly running command in terminal
  actual=t(matrix(c(0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0),ncol=samplesize))
  extract=matrix(unlist(sim$genomes),nrow=samplesize)
  expect_equal(actual,extract)
})


