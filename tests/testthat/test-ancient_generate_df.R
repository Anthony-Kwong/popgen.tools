test_that("ancient_generate_df function works",{
  
  #test 1: random impute ----
  
  #demographic params
  mu=1.5e-8
  recomb_rate=1e-9
  Ne=10000
  genome_length=1e6
  samplesize=200 #try 100 afterwards
  s=0.1
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"
  nsim=5
  ID = seq(6,10,by=1)
  
  #population tree params
  seed = 1
  demes = 2
  sample_dist = c(190,10)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.05
  asc_index = list(c(191,192),c(193,194),c(195,196),c(197,198),c(199,200))
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "random"
  denoise_method = "cluster"
  max_clus = 0.2
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,
                            genome_length=genome_length,samplesize=samplesize,
                            s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep,
                            demes = demes, sample_dist = sample_dist, deme_join = deme_join)
  }
  
  nwins = 5
  
  df<- suppressWarnings(
    ancient_generate_df(l_sim, nwins = nwins,split_type="mut",
                        missing_rate = missing_rate, trans_prop = trans_prop,
                        dmg_rate = dmg_rate, ascertain_indices = asc_index, trim_sim = F, 
                        seed = age_seed,impute_method = impute_method,ID = ID,
                        denoise_method = denoise_method, max_clus = max_clus)
  )

  
  # ancient_sum_stats(sim = l_sim[[1]], nwins = nwins,split_type="mut",
  #                   missing_rate = missing_rate, trans_prop = trans_prop,
  #                   dmg_rate = dmg_rate, ascertain_indices = asc_index, trim_sim = F,
  #                   seed = age_seed,ID = 2)
  
  
  set.seed(age_seed)
  seeds <- sample(.Machine$integer.max, nsim)
  suppressWarnings(
    for(sim_number in 1:nsim){
      batch_ans<-df[sim_number,] 
      batch_ans<-as.numeric(batch_ans)
      single_ans<-ancient_sum_stats(l_sim[[sim_number]],split_type="mut",nwins = nwins,
                                    ID=ID[sim_number], missing_rate = missing_rate,trans_prop = trans_prop,
                                    dmg_rate = dmg_rate,ascertain_indices = asc_index,seed = seeds[sim_number],
                                    impute_method = impute_method, denoise_method = denoise_method,
                                    max_clus = max_clus)
      #single_ans$sweep<-as.factor(single_ans$sweep)
      single_ans<-as.numeric(single_ans)
      expect_equal(batch_ans,single_ans)
    }
  )
  
  #test 2: zero impute ----

  #demographic params
  mu=1.5e-8
  recomb_rate=1e-9
  Ne=10000
  genome_length=1e6
  samplesize=200 #try 100 afterwards
  s=0.1
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"
  nsim=5
  
  #population tree params
  seed = 1
  demes = 2
  sample_dist = c(190,10)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.05
  asc_index = list(c(191,192),c(193,194),c(195,196),c(197,198),c(199,200))
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "zero"
  denoise_method = "fixed_cluster"
  fixed_clus = 10
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,
                            genome_length=genome_length,samplesize=samplesize,
                            s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep,
                            demes = demes, sample_dist = sample_dist, deme_join = deme_join)
  }
  
  nwins = 5
  
  df<- suppressWarnings(
    ancient_generate_df(l_sim, nwins = nwins,split_type="mut",
                        missing_rate = missing_rate, trans_prop = trans_prop,
                        dmg_rate = dmg_rate, ascertain_indices = asc_index, trim_sim = F, 
                        seed = age_seed,impute_method = impute_method,
                        denoise_method = denoise_method, fixed_clus = fixed_clus)
  )
  
  
  # ancient_sum_stats(sim = l_sim[[1]], nwins = nwins,split_type="mut",
  #                   missing_rate = missing_rate, trans_prop = trans_prop,
  #                   dmg_rate = dmg_rate, ascertain_indices = asc_index, trim_sim = F,
  #                   seed = age_seed,ID = 2)
  
  
  set.seed(age_seed)
  seeds <- sample(.Machine$integer.max, nsim)
  suppressWarnings(
    for(sim_number in 1:nsim){
      batch_ans<-df[sim_number,] 
      batch_ans<-as.numeric(batch_ans)
      single_ans<-ancient_sum_stats(l_sim[[sim_number]],split_type="mut",nwins = nwins,
                                    ID=sim_number, missing_rate = missing_rate,trans_prop = trans_prop,
                                    dmg_rate = dmg_rate,ascertain_indices = asc_index,seed = seeds[sim_number],
                                    impute_method = impute_method,
                                    denoise_method = denoise_method, fixed_clus = fixed_clus)
      #single_ans$sweep<-as.factor(single_ans$sweep)
      single_ans<-as.numeric(single_ans)
      expect_equal(batch_ans,single_ans)
    }
  )
  
  #test 3, random_impute , cluster
  
  #demographic params
  mu=1.5e-8
  recomb_rate=1e-9
  Ne=10000
  genome_length=1e6
  samplesize=200 #try 100 afterwards
  s=0.1
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep="hard"
  nsim=5
  ID = seq(6,10,by=1)
  
  #population tree params
  seed = 1
  demes = 2
  sample_dist = c(190,10)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.05
  asc_index = list(c(191,192),c(193,194),c(195,196),c(197,198),c(199,200))
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 5
  impute_method = "random"
  denoise_method = "majority_flip"
  
  l_sim<-list()
  
  for(i in 1:nsim){
    l_sim[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,
                            genome_length=genome_length,samplesize=samplesize,
                            s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep,
                            demes = demes, sample_dist = sample_dist, deme_join = deme_join)
  }
  
  nwins = 5
  
  df<- suppressWarnings(
    ancient_generate_df(l_sim, nwins = nwins,split_type="mut",
                        missing_rate = missing_rate, trans_prop = trans_prop,
                        dmg_rate = dmg_rate, ascertain_indices = asc_index, trim_sim = F, 
                        seed = age_seed,impute_method = impute_method,ID = ID,
                        denoise_method = denoise_method)
  )
  
  
  # ancient_sum_stats(sim = l_sim[[1]], nwins = nwins,split_type="mut",
  #                   missing_rate = missing_rate, trans_prop = trans_prop,
  #                   dmg_rate = dmg_rate, ascertain_indices = asc_index, trim_sim = F,
  #                   seed = age_seed,ID = 2)
  
  
  set.seed(age_seed)
  seeds <- sample(.Machine$integer.max, nsim)
  suppressWarnings(
    for(sim_number in 1:nsim){
      batch_ans<-df[sim_number,] 
      batch_ans<-as.numeric(batch_ans)
      single_ans<-ancient_sum_stats(l_sim[[sim_number]],split_type="mut",nwins = nwins,
                                    ID=ID[sim_number], missing_rate = missing_rate,trans_prop = trans_prop,
                                    dmg_rate = dmg_rate,ascertain_indices = asc_index,seed = seeds[sim_number],
                                    impute_method = impute_method, denoise_method = denoise_method)
      #single_ans$sweep<-as.factor(single_ans$sweep)
      single_ans<-as.numeric(single_ans)
      expect_equal(batch_ans,single_ans)
    }
  )
  
})