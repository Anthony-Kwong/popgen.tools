test_that("Discoal command entered correctly", {

  #input parameters
  mu=1e-8
  recomb_rate=1e-9
  Ne=1000000
  genome_length=1e5
  samplesize=20
  s=0.1
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
  cmd<-sim$cmd

  alpha=no_scientific(2*Ne*s) #scaled strength of selection
  theta=no_scientific(4*Ne*mu*genome_length) #scaled mutation rate
  rho=no_scientific(4*Ne*recomb_rate*genome_length)  #recomb_rate is the probability of a cross over per basepair of sequence being modelled.
  tau= no_scientific(fix/(4*Ne)) #scaled time for fixation

  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(200000),"-t",theta,"-r", rho,
                 "-a", alpha, "-ws", tau)
  expect_equal(cmd,test_cmd)

  seeds=c(1,2)
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  cmd<-sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(200000),"-t",theta,"-r", rho,"-d", seeds[1],seeds[2],
                 "-a", alpha, "-ws", tau)
  expect_equal(cmd,test_cmd)
  
  #test hard sweep command
  sweep="hard"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,fix_generation=fix,sweep=sweep,s=0.2)
  input_cmd=sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(200000),"-t",theta,"-r", rho,
                 "a", alpha, "-ws", tau)
  
  #test soft sweep
  sweep="soft"
  start_freq=0.2
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,fix_generation=fix,sweep=sweep,start_freq=start_freq,s=0.1)
  input_cmd=sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(200000),"-t",theta,"-r", rho,
                 "a", alpha, "-ws", tau, "-f", start_freq)

  #testing the neutral sweep command

  sweep="neutral_fixation"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
  input_cmd=sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(200000),"-t",theta,"-r", rho,
                 "-wn", tau)
  expect_equal(input_cmd,test_cmd)
  
  sweep="neutral"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,sweep=sweep)
  input_cmd=sim$cmd
  test_cmd=paste(discoal_path, no_scientific(samplesize),1,no_scientific(200000),"-t",theta,"-r", rho)
  expect_equal(input_cmd,test_cmd)
})

test_that("Seed extraction successful",{

  #input parameters
  set.seed(123)
  mu=runif(1,1e-8,1e-7)
  recomb_rate=runif(1,1e-10,1e-9)
  Ne=runif(1,1000,1000000)
  genome_length=1e5
  samplesize=runif(1,1,50)
  s=runif(1,1e-3,0.1)
  fix=runif(1,1,5)
  discoal_path="~/work/programs/discoal/discoal"
  seed1=floor(runif(1,1,9e6))
  seed2=floor(runif(1,1,9e6))
  seeds=c(seed1,seed2)
  
  sweep="neutral"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,seed=seeds,sweep=sweep)
  used_seed<-sim$seeds
  expect_equal(used_seed,seeds)
  
  sweep="neutral_fixation"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  used_seed<-sim$seeds
  expect_equal(used_seed,seeds)
  
  sweep="hard"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  used_seed<-sim$seeds
  expect_equal(used_seed,seeds)
  
  sweep="soft"
  start_freq=0.2
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep,start_freq = start_freq)
  used_seed<-sim$seeds
  expect_equal(used_seed,seeds)
})

test_that("Segsites extraction successful",{
  
  #testing for hard sweeps
  
  #input parameters
  mu=2e-8
  recomb_rate=1e-9
  Ne=100
  genome_length=1e6
  samplesize=30
  s=0.1
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  seeds=c(2019,1688)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  extracted_seg=sim$num_seg
  #print(extracted_seg)
  #Ran in terminal to get this
  actual_seg=10
  expect_equal(extracted_seg,actual_seg)
  
  #test case of having no segregating sites
  samplesize=1
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  expect_equal(sim,"Simulation produced no segregating sites")

  #test for neutral simulations
  sweep="neutral"
  mu=2e-6
  recomb_rate=1e-9
  Ne=100
  genome_length=1e5
  samplesize=30
  
  
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,seed=seeds,sweep=sweep)
  extracted_seg=sim$num_seg
  actual_seg=346
  expect_equal(extracted_seg,actual_seg)
  
  #test for neutral fixation
  sweep="neutral_fixation"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,seed=seeds,fix_generation=fix,sweep=sweep)
  extracted_seg=sim$num_seg
  actual_seg=277
  expect_equal(extracted_seg,actual_seg)
  
  freq=0.2
  sweep="soft"
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,s=0.2,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,seed=seeds,fix_generation=fix,sweep=sweep,start_freq = freq)
  extracted_seg=sim$num_seg
  actual_seg=118
  expect_equal(extracted_seg,actual_seg)
})

test_that("Postions extraction successful",{
  #input parameters
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000000
  genome_length=1e5
  samplesize=3
  s=0.1
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  seeds=c(1453,1688)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)
  extracted_pos=sim$pos
  #Ran the command in terminal to get these numbers
  actual_pos=c(0.575357, 0.592501, 0.677530, 0.991324 )
  expect_equal(extracted_pos,actual_pos)
})

test_that("Genome matrix extraction successful",{
  #fix. Code acts funny when there is only one segsite
  
  #input parameters
  mu=2e-7
  recomb_rate=1e-8
  Ne=500
  genome_length=1e5
  samplesize=5
  s=0.1
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  seeds=c(1688,1707)
  sweep="hard"

  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep)

  #Obtained via directly running command in terminal
  actual=cbind(c(1,0,0,0,0),c(0,0,0,1,0),c(0,0,0,1,0),c(1,0,0,0,0),c(0,1,0,0,0))
  extract=matrix(unlist(sim$genomes),nrow=samplesize)
  expect_equal(actual,extract)
  
  
  #input parameters
  mu=2e-7
  recomb_rate=1e-8
  Ne=100
  genome_length=1e5
  samplesize=5
  s=0.2
  fix=2
  seeds=c(1688,1707)
  sweep="soft"
  freq=0.4
  sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,seed=seeds,sweep=sweep,start_freq = freq)
  actual=cbind(c(1,0,0,0,0),c(0,0,1,0,0))
  extract=matrix(unlist(sim$genomes),nrow=samplesize)
  expect_equal(actual,extract)
})


