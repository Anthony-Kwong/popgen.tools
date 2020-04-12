# Downsample test ----

test_that("sum_stats downsampling works",{
  mu=2e-8
  recomb_rate= 0
  Ne=1000
  nBases=1e6
  samplesize=200
  s=0.01
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=5
  id=1
  seed=c(9,8)
  snp_inc = 250
  ds_seed = 52
  ds_prop = 0.1
  
  sim = discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,
                    samplesize=samplesize,s=s,discoal_path=discoal_path,
                    fix_time=fix,sweep=sweep_type,seed=seed)
  
  output = sum_stats(sim = sim,split_type="base",nwins = nwins, ID=id, snp = snp_inc,
                     LD_downsample = T, ds_prop = ds_prop, ds_seed = ds_seed)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  
  snp_dist = abs(sim$pos-0.5)
  center = which.min(snp_dist)
  G = window_trim( sim$genomes ,cen = center, k=floor(snp_inc/2))
  pos_vec<-vector_trim(sim$pos, cen=center, k=floor(snp_inc/2)) 
  G_wins = winsplit_base(G, pos_vec, nwins)

  #check SFS stats
  #check SS were computed correctly
  
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
  h_act<-lapply(G_wins,h_stats) 
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
  
  #check LD stats (downsampled)
  
  set.seed(ds_seed)
  seeds = sample(.Machine$integer.max,nwins)
  down_win_list=purrr::pmap(list(G_wins, p=ds_prop, seed=seeds),
                            downsample_mat)
  
  
  LD_act <- lapply(down_win_list, LD_calc)
  LD_act <- do.call(rbind,LD_act)
  
  LD_avg <- output %>% dplyr::select(LD_avg_1:LD_avg_5) %>% as.numeric()
  LD_max <- output %>% dplyr::select(LD_max_1:LD_max_5) %>% as.numeric()
  w_max <- output %>% dplyr::select(w_max_1:w_max_5) %>% as.numeric()
  Zns <- output %>% dplyr::select(Zns_1:Zns_5) %>% as.numeric()
  
  LD_output= tibble::tibble(LD_avg, LD_max , w_max , Zns)
  expect_equal(LD_act, LD_output)
})

# Hard sweep test without window_trim ----

test_that("sum_stats works",{
  
  #hard sweep test, no transformation. num_seg<snp included
  mu=1.5e-8
  recomb_rate=1e-8
  Ne=1000
  nBases=1e6
  samplesize=200
  s=0.1
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=10
  id=1
  seed=c(1,1)
  
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep_type,seed=seed)
  output<-sum_stats(sim=temp,split_type="mut",nwins = nwins, ID=id, snp=200)
  G<-temp$genomes
  G_wins<-sub_win(G,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  
  
  #check SS were computed correctly
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_10)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_10)
  D_output<-as.numeric(D_output)
  
  LD_act <- lapply(G_wins, LD_calc)
  LD_act <- do.call(rbind,LD_act)
  LD_avg_output <- output %>% dplyr::select(LD_avg_1:LD_avg_10) %>% as.numeric()
  LD_max_output <- output %>% dplyr::select(LD_max_1:LD_max_10) %>% as.numeric()
  w_max_output <- output %>% dplyr::select(w_max_1:w_max_10) %>% as.numeric()
  Zns_output <- output %>% dplyr::select(Zns_1:Zns_10) %>% as.numeric()
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }

  h1_output = output %>% dplyr::select(h1_1:h1_10)
  h2_output = output %>% dplyr::select(h2_1:h2_10)
  h12_output = output %>% dplyr::select(h12_1:h12_10)
  h123_output = output %>% dplyr::select(h123_1:h123_10)
  
  expect_equal(H_act,H_output)
  expect_equal(D_act,D_output)
  expect_equal(LD_act$LD_avg,LD_avg_output)
  expect_equal(LD_act$LD_max,LD_max_output)
  expect_equal(LD_act$w_max,w_max_output)
  expect_equal(LD_act$Zns,Zns_output)
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
})

## Hard simulation test with window_trim ----

test_that("sum_stats works",{
  
  #hard sweep test, no transformation, num_seg>snp included (testing trimming)
  mu=1.5e-8
  recomb_rate=1e-8
  Ne=1000
  nBases=1e6
  samplesize=200
  s=0.1
  fix=1
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=10
  id=1
  seed=c(1,2)
  snp_inc=200
  
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep_type,seed=seed)
  output<-sum_stats(sim=temp,split_type="mut",nwins = nwins, ID=id, snp=snp_inc)
  
  #trim matrix
  snp_dist<-abs(temp$pos-0.5)
  center<-which.min(snp_dist)
  G<-temp$genomes
  G<-window_trim(temp$genomes,cen=center,k=floor(snp_inc/2))
  G_wins<-sub_win(G,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)
  
  #check SS were computed correctly
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_10)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_10)
  D_output<-as.numeric(D_output)
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_10)
  h2_output = output %>% dplyr::select(h2_1:h2_10)
  h12_output = output %>% dplyr::select(h12_1:h12_10)
  h123_output = output %>% dplyr::select(h123_1:h123_10)
  
  expect_equal(H_act,H_output)
  expect_equal(D_act,D_output)
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  #LD stats 
  
  LD_act <- lapply(G_wins, LD_calc)
  LD_act <- do.call(rbind,LD_act)
  LD_avg_output <- output %>% dplyr::select(LD_avg_1:LD_avg_10) %>% as.numeric()
  LD_max_output <- output %>% dplyr::select(LD_max_1:LD_max_10) %>% as.numeric()
  w_max_output <- output %>% dplyr::select(w_max_1:w_max_10) %>% as.numeric()
  Zns_output <- output %>% dplyr::select(Zns_1:Zns_10) %>% as.numeric()
  
  expect_equal(LD_act$LD_avg,LD_avg_output)
  expect_equal(LD_act$LD_max,LD_max_output)
  expect_equal(LD_act$w_max,w_max_output)
  expect_equal(LD_act$Zns,Zns_output)
  
})

## Norm transformation across window test ---- 


test_that("sum_stat works",{
  #neutral test, norm transformation
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000
  nBases=1e6
  samplesize=250
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="neutral"
  nwins=10
  id=12
  seed=c(13,82)
  snp_inc=450

  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,discoal_path=discoal_path,sweep=sweep_type,seed=seed)
  output<-sum_stats(sim=temp,split_type="mut",nwins = nwins, ID=id, snp=snp_inc,form="wide",fun="norm")
  
  #trim matrix
  snp_dist<-abs(temp$pos-0.5)
  center<-which.min(snp_dist)
  G<-temp$genomes
  G<-window_trim(temp$genomes,cen=center,k=floor(snp_inc/2))
  G_wins<-sub_win(G,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)

  #check SS were computed correctly
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_act <- norm_vec(H_act)
  H_output<-output %>% dplyr::select(H_1:H_10) 
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_act<-norm_vec(D_act)
  D_output<-output %>% dplyr::select(D_1:D_10)
  D_output<-as.numeric(D_output)
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_10)
  h2_output = output %>% dplyr::select(h2_1:h2_10)
  h12_output = output %>% dplyr::select(h12_1:h12_10)
  h123_output = output %>% dplyr::select(h123_1:h123_10)
  
  expect_equal(H_act,H_output)
  expect_equal(D_act,D_output)
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  #LD stats
  
  LD_act <- lapply(G_wins, LD_calc)
  LD_act <- do.call(rbind,LD_act)
  LD_avg_act <- LD_act$LD_avg %>% norm_vec()
  LD_max_act <- LD_act$LD_max %>% norm_vec()
  w_max_act <- LD_act$w_max %>% norm_vec()
  Zns_act <- LD_act$Zns %>% norm_vec()
  
  LD_avg_output <- output %>% dplyr::select(LD_avg_1:LD_avg_10) %>% as.numeric() %>% norm_vec()
  LD_max_output <- output %>% dplyr::select(LD_max_1:LD_max_10) %>% as.numeric() %>% norm_vec()
  w_max_output <- output %>% dplyr::select(w_max_1:w_max_10) %>% as.numeric() %>% norm_vec()
  Zns_output <- output %>% dplyr::select(Zns_1:Zns_10) %>% as.numeric() %>% norm_vec()
  
  
  expect_equal(LD_avg_act,LD_avg_output)
  expect_equal(LD_max_act,LD_max_output)
  expect_equal(w_max_act,w_max_output)
  expect_equal(Zns_act,Zns_output)
})

## Split windows by base test ---- 

test_that("sum_stats works",{
  
  #neutral test, split by base
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000
  nBases=1e6
  samplesize=250
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="neutral"
  nwins=10
  id=2
  seed=c(14,81)
  snp_inc=380
  split_type="base"

  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,discoal_path=discoal_path,sweep=sweep_type,seed=seed)
  output<-sum_stats(sim=temp,split_type=split_type,nwins = nwins, ID=id, snp=snp_inc,form="wide",fun="none")

  #trim matrix
  snp_dist<-abs(temp$pos-0.5)
  center<-which.min(snp_dist)
  G<-window_trim(temp$genomes,cen=center,k=floor(snp_inc/2))
  positions = vector_trim(temp$pos,center,k=floor(snp_inc/2))
  G_wins<-winsplit_base(G,pos=positions,nwins)

  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(0,output$s_coef)
  expect_equal(0,output$bottle_time1)
  expect_equal(0,output$bottle_time2)
  expect_equal(1,output$bottle_size1)
  expect_equal(1,output$bottle_size2)

  #check SS were computed correctly

  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_10)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_10)
  D_output<-as.numeric(D_output)
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_10)
  h2_output = output %>% dplyr::select(h2_1:h2_10)
  h12_output = output %>% dplyr::select(h12_1:h12_10)
  h123_output = output %>% dplyr::select(h123_1:h123_10)
  
  expect_equal(H_act,H_output)
  expect_equal(D_act,D_output)
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  #LD stats 
  
  LD_act <- lapply(G_wins, LD_calc)
  LD_act <- do.call(rbind,LD_act)
  LD_avg_output <- output %>% dplyr::select(LD_avg_1:LD_avg_10) %>% as.numeric()
  LD_max_output <- output %>% dplyr::select(LD_max_1:LD_max_10) %>% as.numeric()
  w_max_output <- output %>% dplyr::select(w_max_1:w_max_10) %>% as.numeric()
  
  expect_equal(LD_act$LD_avg,LD_avg_output)
  expect_equal(LD_act$LD_max,LD_max_output)
  expect_equal(LD_act$w_max,w_max_output)

})

## Bottleneck tests----

test_that("sum_stats works",{
  
  #neutral test, split by base
  mu=2e-8
  recomb_rate=1e-8
  Ne=1000
  nBases=1e6
  samplesize=250
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="neutral"
  nwins=10
  id=2
  seed=c(15,101)
  snp_inc=420
  split_type="base"
  bottleneck=tibble::tibble(size=c(0.2,1),time=c(5,10))
  
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,
                    genome_length=nBases,samplesize=samplesize,
                    discoal_path=discoal_path,sweep=sweep_type,
                    seed=seed, popsize_changes = bottleneck)
  output<-sum_stats(sim=temp,split_type=split_type,nwins = nwins, ID=id, snp=snp_inc,form="wide",fun="none")
  
  #trim matrix
  snp_dist<-abs(temp$pos-0.5)
  center<-which.min(snp_dist)
  G<-window_trim(temp$genomes,cen=center,k=floor(snp_inc/2))
  positions = vector_trim(temp$pos,center,k=floor(snp_inc/2))
  G_wins<-winsplit_base(G,pos=positions,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(0,output$s_coef)
  expect_equal(bottleneck$time[1],output$bottle_time1)
  expect_equal(bottleneck$time[2],output$bottle_time2)
  expect_equal(bottleneck$size[1],output$bottle_size1)
  expect_equal(bottleneck$size[2],output$bottle_size2)
  
  #check SS were computed correctly
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output %>% dplyr::select(H_1:H_10)
  H_output<-as.numeric(H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output %>% dplyr::select(D_1:D_10)
  D_output<-as.numeric(D_output)
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  h1_output = output %>% dplyr::select(h1_1:h1_10)
  h2_output = output %>% dplyr::select(h2_1:h2_10)
  h12_output = output %>% dplyr::select(h12_1:h12_10)
  h123_output = output %>% dplyr::select(h123_1:h123_10)
  
  expect_equal(H_act,H_output)
  expect_equal(D_act,D_output)
  expect_equal(as.numeric(h1_output),haplo_act[,1])
  expect_equal(as.numeric(h2_output),haplo_act[,2])
  expect_equal(as.numeric(h12_output),haplo_act[,3])
  expect_equal(as.numeric(h123_output),haplo_act[,4])
  
  #LD stats 
  
  LD_act <- lapply(G_wins, LD_calc)
  LD_act <- do.call(rbind,LD_act)
  LD_avg_output <- output %>% dplyr::select(LD_avg_1:LD_avg_10) %>% as.numeric()
  LD_max_output <- output %>% dplyr::select(LD_max_1:LD_max_10) %>% as.numeric()
  w_max_output <- output %>% dplyr::select(w_max_1:w_max_10) %>% as.numeric()
  
  expect_equal(LD_act$LD_avg,LD_avg_output)
  expect_equal(LD_act$LD_max,LD_max_output)
  expect_equal(LD_act$w_max,w_max_output)
  
})

