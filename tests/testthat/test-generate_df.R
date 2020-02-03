test_that("generate_df works",{
  
  #hard sweep simulation test, checking fun="norm"
  
  #input parameters
  mu=1.5e-8
  recomb_rate=1e-9
  Ne=1000
  genome_length=1e6
  samplesize=10
  s=0.01
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"
  nsim=5
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
  }
  
  snp_include=50
  nwins=3
  
  df<-generate_df(l_sim,win_split = nwins,snp = snp_include,form="wide",fun="norm")
  
  #check row 3 numbers
  
  for(index in 1:nsim){
    batch_ans<-df[index,] 
    batch_ans<-as.numeric(batch_ans)
    single_ans<-sum_stats(l_sim[[index]],win_split = nwins, ID=index, snp=snp_include,fun="norm")
    single_ans$sweep<-as.factor(single_ans$sweep)
    single_ans<-as.numeric(single_ans)
    expect_equal(batch_ans,single_ans)
  }
  
## neutral test with default params for form and fun---

  #input parameters
  mu=1e-8
  recomb_rate=1e-9
  Ne=10000
  genome_length=1e5
  samplesize=20
  discoal_path="~/work/programs/discoal/discoal"
  sweep="neutral"
  nsim=5
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,discoal_path=discoal_path,sweep=sweep)
  }
  
  snp_include=100
  nwins=4
  
  df<-generate_df(l_sim,win_split = nwins,snp = snp_include)
  
  for(index in 1:nsim){
    batch_ans<-df[index,] 
    batch_ans<-as.numeric(batch_ans)
    single_ans<-sum_stats(l_sim[[index]],win_split = nwins, ID=index, snp=snp_include)
    single_ans$sweep<-as.factor(single_ans$sweep)
    single_ans<-as.numeric(single_ans)
    expect_equal(batch_ans,single_ans)
  }
  
  ## hard sweep test with default params for form and fun
  
  #hard sweep simulation test, checking fun="norm"
  
  #input parameters
  mu=1.6e-8
  recomb_rate=1e-9
  Ne=1000
  genome_length=1e6
  samplesize=10
  s=0.005
  fix=2
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"
  nsim=10
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
  }
  
  snp_include=50
  nwins=5
  
  df<-generate_df(l_sim,win_split = nwins,snp = snp_include,form="wide")
  
  #check row 3 numbers
  
  for(index in 1:nsim){
    batch_ans<-df[index,] 
    batch_ans<-as.numeric(batch_ans)
    single_ans<-sum_stats(l_sim[[index]],win_split = nwins, ID=index, snp=snp_include)
    single_ans$sweep<-as.factor(single_ans$sweep)
    single_ans<-as.numeric(single_ans)
    expect_equal(batch_ans,single_ans)
  }
  
})

