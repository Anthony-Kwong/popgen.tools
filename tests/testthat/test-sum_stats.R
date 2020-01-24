test_that("sum_stats works",{
  
  #convert sweep to neutral if the selection coefficient is 0. 
  
  #hard sweep test
  mu=1.5e-8
  recomb_rate=1e-8
  Ne=1000
  nBases=1e6
  samplesize=200
  s=0.1
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep_type)
  sum_stats(sim=temp,win_split = 10, ID=1, snp=200)
})

