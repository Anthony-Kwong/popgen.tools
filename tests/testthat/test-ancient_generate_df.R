test_that("ancient_generate_df function works",{
  
  #demographic params
  mu=1.5e-8
  recomb_rate=1e-9
  Ne=1000
  genome_length=1e6
  samplesize=25
  s=0.01
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"
  nsim=5
  
  #population tree params
  seed = 1
  demes = 2
  sample_dist = c(23,2)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.05
  asc_index = c(23,24,25)
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,
                            genome_length=genome_length,samplesize=samplesize,
                            s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep,
                            demes = demes, sample_dist = sample_dist, deme_join = deme_join)
  }
  
  nwins = 5
  
  df<-ancient_generate_df(l_sim, nwins = nwins,split_type="mut",
                          missing_rate = missing_rate, trans_prop = trans_prop,
                          dmg_rate = dmg_rate, index = asc_index, trim_sim = F, 
                          seed = age_seed)
  
  set.seed(age_seed)
  seeds <- sample(.Machine$integer.max, nsim)
  suppressWarnings(
    for(sim_number in 1:nsim){
      batch_ans<-df[sim_number,] 
      batch_ans<-as.numeric(batch_ans)
      single_ans<-ancient_sum_stats(l_sim[[sim_number]],split_type="mut",nwins = nwins,
                                    ID=sim_number, missing_rate = missing_rate,trans_prop = trans_prop,
                                    dmg_rate = dmg_rate,index = asc_index,seed = seeds[sim_number])
      #single_ans$sweep<-as.factor(single_ans$sweep)
      single_ans<-as.numeric(single_ans)
      expect_equal(batch_ans,single_ans)
    }
  )
  
})