test_that("ancient_sum_stats works",{
  
  # test 1: random impute ----
  
  #simulation params
  mu=2e-8
  recomb_rate= 0
  Ne=10000
  nBases=1e6
  samplesize=200
  s=0.001
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=5
  id=1
  seed=c(9,8)
  snp_inc = 250
  demes = 2
  sample_dist = c(196,4)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.05
  asc_index = list(c(195,196),c(197,198),c(199,200))
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "random"
  
  
  sim = discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,
                    samplesize=samplesize,s=s,discoal_path=discoal_path,
                    fix_time=fix,sweep=sweep_type,seed=seed,
                    demes = demes,sample_dist = sample_dist,deme_join = deme_join)
  
  output = suppressWarnings(
    ancient_sum_stats(sim = sim,split_type="mut",nwins = nwins, ID=id, snp = snp_inc,
                      ascertain_indices = asc_index, missing_rate = missing_rate,trans_prop = trans_prop
                      ,seed = age_seed,impute_method = impute_method)
  )
  
  #check non SS information
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  expect_equal(output$impute_method,impute_method)
  
  
  #aging DNA
  raw_G <- sim$genomes
  
  num_pairs <- length(asc_index)
  asc_sites_list <- lapply(seq(num_pairs), function(d){ascertain_bias(raw_G, asc_index[[d]])[[2]] } )
  asc_sites <- asc_sites_list %>%
    unlist() %>% 
    unique() %>%
    sort()
  asc_rows <- asc_index %>%
    unlist()
  asc_G = raw_G[-c(asc_rows),asc_sites]
  
  dmg_G <- age_DNA(G = asc_G, missing_rate = missing_rate,
                   trans_prop = trans_prop,dmg_rate = dmg_rate, seed = age_seed)
  dmg_G <- rm_nonpoly_cols(dmg_G)
  
  set.seed(age_seed)
  imp_G <- random_impute(dmg_G[[1]])
  imp_pos <- sim$pos[asc_sites]
  
  #remove non-polymorphic sites
  rm_G = rm_nonpoly_cols(imp_G)
  final_G = rm_G[[1]]
  final_pos = imp_pos[c(rm_G[[2]])]
  
  #split windows based by ~equal SNPs
  split_wins = sub_win(final_G,nwins)
  win_list= split_wins$windows
  
  #compute the length of each block in bases
  j = seq(2, nwins + 1)
  base_lengths = sapply(j, function(j){ final_pos[ split_wins$bounds[j] ] - 
      final_pos[ split_wins$bounds[j-1] ]})
  base_lengths = as.data.frame(base_lengths) %>% t() %>% as.numeric()
  base_output = output %>% dplyr::select(block_base_length_1:block_base_length_5) %>% as.numeric()
  expect_equal(base_output, base_lengths)
  
  #compute number of SNPs in each block. This is more for checking purposes.
  snp_output = output %>% dplyr::select(block_snp_length_1:block_snp_length_5) %>% as.numeric()
  snp_act = sapply(win_list, ncol)
  expect_equal(snp_output, snp_act)
  
  #check SS computation
  
  win_splits = sub_win(final_G,nwins)
  G_wins = win_splits[[1]]
  
  #check SFS stats
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_5)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_5)
  D_output<-as.numeric(D_output)
  
  expect_equal(H_output, H_act)
  expect_equal(D_output, D_act)
  
  #check haplotype stats
  set.seed(age_seed)
  pseudo_G = pseudo_hap(final_G, seed = age_seed)
  
  G_wins = sub_win(pseudo_G,nwins)
  hap_win_list= G_wins$windows
  
  h_act<-lapply(hap_win_list,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_5)
  h2_output = output %>% dplyr::select(h2_1:h2_5)
  h12_output = output %>% dplyr::select(h12_1:h12_5)
  h123_output = output %>% dplyr::select(h123_1:h123_5)
  
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  # test 2: zero impute ----
  
  #simulation params
  mu=2e-8
  recomb_rate= 0
  Ne=10000
  nBases=1e6
  samplesize=200
  s=0.001
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=5
  id=1
  seed=c(9,8)
  snp_inc = 250
  demes = 2
  sample_dist = c(196,4)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.05
  asc_index = list(c(195,196),c(197,198),c(199,200))
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "zero"
  
  
  sim = discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,
                    samplesize=samplesize,s=s,discoal_path=discoal_path,
                    fix_time=fix,sweep=sweep_type,seed=seed,
                    demes = demes,sample_dist = sample_dist,deme_join = deme_join)
  
  output = suppressWarnings(
    ancient_sum_stats(sim = sim,split_type="mut",nwins = nwins, ID=id, snp = snp_inc,
                      ascertain_indices = asc_index, missing_rate = missing_rate,trans_prop = trans_prop
                      ,seed = age_seed,impute_method = impute_method)
  )
  
  #check non SS information
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  expect_equal(output$impute_method,impute_method)
  
  
  #aging DNA
  raw_G <- sim$genomes
  
  num_pairs <- length(asc_index)
  asc_sites_list <- lapply(seq(num_pairs), function(d){ascertain_bias(raw_G, asc_index[[d]])[[2]] } )
  asc_sites <- asc_sites_list %>%
    unlist() %>% 
    unique() %>%
    sort()
  asc_rows <- asc_index %>%
    unlist()
  asc_G = raw_G[-c(asc_rows),asc_sites]
  
  dmg_G <- age_DNA(G = asc_G, missing_rate = missing_rate,
                   trans_prop = trans_prop,dmg_rate = dmg_rate, seed = age_seed)
  dmg_G <- rm_nonpoly_cols(dmg_G)
  
  set.seed(age_seed)
  imp_G <- zero_impute(dmg_G[[1]])
  imp_pos <- sim$pos[asc_sites]
  
  #remove non-polymorphic sites
  rm_G = rm_nonpoly_cols(imp_G)
  final_G = rm_G[[1]]
  final_pos = imp_pos[c(rm_G[[2]])]
  
  #split windows based by ~equal SNPs
  split_wins = sub_win(final_G,nwins)
  win_list= split_wins$windows
  
  #compute the length of each block in bases
  j = seq(2, nwins + 1)
  base_lengths = sapply(j, function(j){ final_pos[ split_wins$bounds[j] ] - 
      final_pos[ split_wins$bounds[j-1] ]})
  base_lengths = as.data.frame(base_lengths) %>% t() %>% as.numeric()
  base_output = output %>% dplyr::select(block_base_length_1:block_base_length_5) %>% as.numeric()
  expect_equal(base_output, base_lengths)
  
  #compute number of SNPs in each block. This is more for checking purposes.
  snp_output = output %>% dplyr::select(block_snp_length_1:block_snp_length_5) %>% as.numeric()
  snp_act = sapply(win_list, ncol)
  expect_equal(snp_output, snp_act)
  
  #check SS computation
  
  win_splits = sub_win(final_G,nwins)
  G_wins = win_splits[[1]]
  
  #check SFS stats
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_5)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_5)
  D_output<-as.numeric(D_output)
  
  expect_equal(H_output, H_act)
  expect_equal(D_output, D_act)
  
  #check haplotype stats
  set.seed(age_seed)
  pseudo_G = pseudo_hap(final_G, seed = age_seed)
  
  G_wins = sub_win(pseudo_G,nwins)
  hap_win_list= G_wins$windows
  
  h_act<-lapply(hap_win_list,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_5)
  h2_output = output %>% dplyr::select(h2_1:h2_5)
  h12_output = output %>% dplyr::select(h12_1:h12_5)
  h123_output = output %>% dplyr::select(h123_1:h123_5)
  
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  # test 3: random impute ----
  
  #simulation params
  mu=1.5e-8
  recomb_rate= 0
  Ne=10000
  nBases=1e6
  samplesize=220
  s=0.1
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=5
  id=1
  seed=c(9,8)
  demes = 2
  sample_dist = c(200,20)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.1
  asc_index = lapply(seq(180,219,by=2), function(d){c(d,d+1)})
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "random"
  
  
  sim = discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,
                    samplesize=samplesize,s=s,discoal_path=discoal_path,
                    fix_time=fix,sweep=sweep_type,seed=seed,
                    demes = demes,sample_dist = sample_dist,deme_join = deme_join)
  
  output = suppressWarnings(
    ancient_sum_stats(sim = sim,split_type="mut",nwins = nwins, ID=id,
                      ascertain_indices = asc_index, missing_rate = missing_rate,trans_prop = trans_prop
                      ,seed = age_seed,impute_method = impute_method)
  )
  
  #check non SS information
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  expect_equal(output$impute_method,impute_method)
  
  
  #aging DNA
  raw_G <- sim$genomes
  
  num_pairs <- length(asc_index)
  asc_sites_list <- lapply(seq(num_pairs), function(d){ascertain_bias(raw_G, asc_index[[d]])[[2]] } )
  asc_sites <- asc_sites_list %>%
    unlist() %>% 
    unique() %>%
    sort()
  asc_rows <- asc_index %>%
    unlist()
  asc_G = raw_G[-c(asc_rows),asc_sites]
  
  dmg_G <- age_DNA(G = asc_G, missing_rate = missing_rate,
                   trans_prop = trans_prop,dmg_rate = dmg_rate, seed = age_seed)
  dmg_G <- rm_nonpoly_cols(dmg_G)
  
  set.seed(age_seed)
  imp_G <- random_impute(dmg_G[[1]])
  imp_pos <- sim$pos[asc_sites]
  
  #remove non-polymorphic sites
  rm_G = rm_nonpoly_cols(imp_G)
  final_G = rm_G[[1]]
  final_pos = imp_pos[c(rm_G[[2]])]
  
  #split windows based by ~equal SNPs
  split_wins = sub_win(final_G,nwins)
  win_list= split_wins$windows
  
  #compute the length of each block in bases
  j = seq(2, nwins + 1)
  base_lengths = sapply(j, function(j){ final_pos[ split_wins$bounds[j] ] - 
      final_pos[ split_wins$bounds[j-1] ]})
  base_lengths = as.data.frame(base_lengths) %>% t() %>% as.numeric()
  base_output = output %>% dplyr::select(block_base_length_1:block_base_length_5) %>% as.numeric()
  expect_equal(base_output, base_lengths)
  
  #compute number of SNPs in each block. This is more for checking purposes.
  snp_output = output %>% dplyr::select(block_snp_length_1:block_snp_length_5) %>% as.numeric()
  snp_act = sapply(win_list, ncol)
  expect_equal(snp_output, snp_act)
  
  #check SS computation
  
  win_splits = sub_win(final_G,nwins)
  G_wins = win_splits[[1]]
  
  #check SFS stats
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_5)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_5)
  D_output<-as.numeric(D_output)
  
  expect_equal(H_output, H_act)
  expect_equal(D_output, D_act)
  
  #check haplotype stats
  set.seed(age_seed)
  pseudo_G = pseudo_hap(final_G, seed = age_seed)
  
  G_wins = sub_win(pseudo_G,nwins)
  hap_win_list= G_wins$windows
  
  h_act<-lapply(hap_win_list,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_5)
  h2_output = output %>% dplyr::select(h2_1:h2_5)
  h12_output = output %>% dplyr::select(h12_1:h12_5)
  h123_output = output %>% dplyr::select(h123_1:h123_5)
  
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  #test 4, denoise_method = "cluster" ----
  
  #simulation params
  mu=1.5e-8
  recomb_rate= 0
  Ne=10000
  nBases=1e6
  samplesize=220
  s=0.1
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=5
  id=1
  seed=c(5,8)
  demes = 2
  sample_dist = c(200,20)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)
  
  #DNA aging params
  missing_rate = 0.1
  asc_index = lapply(seq(180,219,by=2), function(d){c(d,d+1)})
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "random"
  denoise_method = "cluster"
  max_clus = 0.25
  
  
  sim = discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,
                    samplesize=samplesize,s=s,discoal_path=discoal_path,
                    fix_time=fix,sweep=sweep_type,seed=seed,
                    demes = demes,sample_dist = sample_dist,deme_join = deme_join)
  
  output = suppressWarnings(
    ancient_sum_stats(sim = sim,split_type="mut",nwins = nwins, ID=id,
                      ascertain_indices = asc_index, missing_rate = missing_rate,trans_prop = trans_prop
                      ,seed = age_seed,impute_method = impute_method,
                      denoise_method = denoise_method, max_clus = max_clus)
  )
  
  #check non SS information
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  expect_equal(output$impute_method,impute_method)
  
  
  #aging DNA
  raw_G <- sim$genomes
  
  num_pairs <- length(asc_index)
  asc_sites_list <- lapply(seq(num_pairs), function(d){ascertain_bias(raw_G, asc_index[[d]])[[2]] } )
  asc_sites <- asc_sites_list %>%
    unlist() %>% 
    unique() %>%
    sort()
  asc_rows <- asc_index %>%
    unlist()
  asc_G = raw_G[-c(asc_rows),asc_sites]
  
  dmg_G <- age_DNA(G = asc_G, missing_rate = missing_rate,
                   trans_prop = trans_prop,dmg_rate = dmg_rate, seed = age_seed)
  dmg_G <- rm_nonpoly_cols(dmg_G)
  
  set.seed(age_seed)
  imp_G <- random_impute(dmg_G[[1]])
  imp_pos <- sim$pos[asc_sites]
  
  #remove non-polymorphic sites
  rm_G = rm_nonpoly_cols(imp_G)
  final_G = rm_G[[1]]
  final_pos = imp_pos[c(rm_G[[2]])]
  
  #split windows based by ~equal SNPs
  split_wins = sub_win(final_G,nwins)
  win_list= split_wins$windows
  
  #compute the length of each block in bases
  j = seq(2, nwins + 1)
  base_lengths = sapply(j, function(j){ final_pos[ split_wins$bounds[j] ] - 
      final_pos[ split_wins$bounds[j-1] ]})
  base_lengths = as.data.frame(base_lengths) %>% t() %>% as.numeric()
  base_output = output %>% dplyr::select(block_base_length_1:block_base_length_5) %>% as.numeric()
  expect_equal(base_output, base_lengths)
  
  #compute number of SNPs in each block. This is more for checking purposes.
  snp_output = output %>% dplyr::select(block_snp_length_1:block_snp_length_5) %>% as.numeric()
  snp_act = sapply(win_list, ncol)
  expect_equal(snp_output, snp_act)
  
  #check SS computation
  
  win_splits = sub_win(final_G,nwins)
  G_wins = win_splits[[1]]
  
  #check SFS stats
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_5)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_5)
  D_output<-as.numeric(D_output)
  
  expect_equal(H_output, H_act)
  expect_equal(D_output, D_act)
  
  #check haplotype stats
  set.seed(age_seed)
  pseudo_G = pseudo_hap(final_G, seed = age_seed)
  
  G_wins = sub_win(pseudo_G,nwins)
  hap_win_list= G_wins$windows
  
  win_clus_vec = lapply(hap_win_list, function(G){clus_hap(G, max_clus = round( nrow(G)*max_clus) ) })
  h_values = lapply(win_clus_vec,clus_hstats)
  
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_values[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_5)
  h2_output = output %>% dplyr::select(h2_1:h2_5)
  h12_output = output %>% dplyr::select(h12_1:h12_5)
  h123_output = output %>% dplyr::select(h123_1:h123_5)
  
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  #if haplotyp stats disagree, check the max_clus argument for the clustering
  
  #test 5, denoise_method = "col_flip"
  
  
  #simulation params
  mu=1.5e-8
  recomb_rate= 0
  Ne=10000
  nBases=1e6
  samplesize=220
  s=0.1
  fix=0
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=5
  id=1
  seed=c(9,8)
  demes = 2
  sample_dist = c(200,20)
  deme_join = tibble::tibble(time = 50000, pop1 = 0, pop2 = 1)

  #DNA aging params
  missing_rate = 0.1
  asc_index = lapply(seq(180,219,by=2), function(d){c(d,d+1)})
  trans_prop = 0.77
  dmg_rate = 0.05
  age_seed = 4
  impute_method = "zero"
  denoise_method = "majority_flip"


  sim = discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,
                    samplesize=samplesize,s=s,discoal_path=discoal_path,
                    fix_time=fix,sweep=sweep_type,seed=seed,
                    demes = demes,sample_dist = sample_dist,deme_join = deme_join)

  output = suppressWarnings(
    ancient_sum_stats(sim = sim,split_type="mut",nwins = nwins, ID=id,
                      ascertain_indices = asc_index, missing_rate = missing_rate,trans_prop = trans_prop
                      ,seed = age_seed,impute_method = impute_method, denoise_method = denoise_method)
  )

  #check non SS information
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  expect_equal(output$impute_method,impute_method)


  #aging DNA
  raw_G <- sim$genomes

  num_pairs <- length(asc_index)
  asc_sites_list <- lapply(seq(num_pairs), function(d){ascertain_bias(raw_G, asc_index[[d]])[[2]] } )
  asc_sites <- asc_sites_list %>%
    unlist() %>%
    unique() %>%
    sort()
  asc_rows <- asc_index %>%
    unlist()
  asc_G = raw_G[-c(asc_rows),asc_sites]

  dmg_G <- age_DNA(G = asc_G, missing_rate = missing_rate,
                   trans_prop = trans_prop,dmg_rate = dmg_rate, seed = age_seed)
  dmg_G <- rm_nonpoly_cols(dmg_G)

  set.seed(age_seed)
  imp_G <- zero_impute(dmg_G[[1]])
  imp_pos <- sim$pos[asc_sites]

  #remove non-polymorphic sites
  rm_G = rm_nonpoly_cols(imp_G)
  final_G = rm_G[[1]]
  final_pos = imp_pos[c(rm_G[[2]])]

  #split windows based by ~equal SNPs
  split_wins = sub_win(final_G,nwins)
  win_list= split_wins$windows

  #compute the length of each block in bases
  j = seq(2, nwins + 1)
  base_lengths = sapply(j, function(j){ final_pos[ split_wins$bounds[j] ] -
      final_pos[ split_wins$bounds[j-1] ]})
  base_lengths = as.data.frame(base_lengths) %>% t() %>% as.numeric()
  base_output = output %>% dplyr::select(block_base_length_1:block_base_length_5) %>% as.numeric()
  expect_equal(base_output, base_lengths)

  #compute number of SNPs in each block. This is more for checking purposes.
  snp_output = output %>% dplyr::select(block_snp_length_1:block_snp_length_5) %>% as.numeric()
  snp_act = sapply(win_list, ncol)
  expect_equal(snp_output, snp_act)

  #check SS computation

  win_splits = sub_win(final_G,nwins)
  G_wins = win_splits[[1]]

  #check SFS stats

  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")

  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_5)
  H_output<-as.numeric(H_output)

  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_5)
  D_output<-as.numeric(D_output)

  expect_equal(H_output, H_act)
  expect_equal(D_output, D_act)

  #check haplotype stats
  set.seed(age_seed)
  pseudo_G = pseudo_hap(final_G, seed = age_seed)
  pseudo_G = majority_flip(pseudo_G)

  G_wins = sub_win(pseudo_G,nwins)
  hap_win_list= G_wins$windows

  h_act<-lapply(hap_win_list,h_stats)
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)

  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }

  h1_output = output %>% dplyr::select(h1_1:h1_5)
  h2_output = output %>% dplyr::select(h2_1:h2_5)
  h12_output = output %>% dplyr::select(h12_1:h12_5)
  h123_output = output %>% dplyr::select(h123_1:h123_5)

  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
})