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
  
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep_type,seed=seed)
  output<-sum_stats(sim=temp,win_split = nwins, ID=id, snp=200)
  G<-temp$genomes
  G_wins<-sub_win(G,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  
  #check SS were computed correctly
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output[4:13] 
  H_output<-as.numeric(H_output)
  
  expect_equal(H_act,H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output[14:23]
  D_output<-as.numeric(D_output)
  
  expect_equal(D_act,D_output)
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  expect_equal(as.numeric(output[24:33]),haplo_act[,1])
  expect_equal(as.numeric(output[34:43]),haplo_act[,2])
  expect_equal(as.numeric(output[44:53]),haplo_act[,3])
  expect_equal(as.numeric(output[54:63]),haplo_act[,4])
  
})

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
  
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep_type,seed=seed)
  output<-sum_stats(sim=temp,win_split = nwins, ID=id, snp=snp_inc)
  
  #trim matrix
  snp_dist<-abs(temp$pos-0.5)
  center<-which.min(snp_dist)
  G<-temp$genomes
  G<-window_trim(temp$genomes,cen=center,k=floor(snp_inc/2))
  G_wins<-sub_win(G,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)
  expect_equal(s,output$s_coef)
  
  #check SS were computed correctly
  
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  H_output<-output[4:13] 
  H_output<-as.numeric(H_output)
  
  expect_equal(H_act,H_output)
  
  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act)
  D_output<-output[14:23]
  D_output<-as.numeric(D_output)
  
  expect_equal(D_act,D_output)
  
  h_act<-lapply(G_wins,h_stats) 
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)
  
  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }
  
  expect_equal(as.numeric(output[24:33]),haplo_act[,1])
  expect_equal(as.numeric(output[34:43]),haplo_act[,2])
  expect_equal(as.numeric(output[44:53]),haplo_act[,3])
  expect_equal(as.numeric(output[54:63]),haplo_act[,4])
  
})

test_that("sum_stat works",{
  
  #test if sweep is converted to neutral if s=0 
  mu=2e-8
  recomb_rate=1e-9
  Ne=1000
  nBases=1e6
  samplesize=50
  s=0
  fix=4
  discoal_path="~/work/programs/discoal/discoal"
  sweep_type="hard"
  nwins=10
  id=1
  temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,sweep=sweep_type)
  #This was meant to produce a warning. 
  expect_equal("neutral",temp$sweep)
})

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
  output<-sum_stats(sim=temp,win_split = nwins, ID=id, snp=snp_inc,form="wide",fun="norm")
  
  #trim matrix
  snp_dist<-abs(temp$pos-0.5)
  center<-which.min(snp_dist)
  G<-temp$genomes
  G<-window_trim(temp$genomes,cen=center,k=floor(snp_inc/2))
  G_wins<-sub_win(G,nwins)
  
  expect_equal(sweep_type,output$sweep)
  expect_equal(id,output$ID)

  #check SS were computed correctly

  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(G_wins, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  #normalise across the values of H and compare with output (already normed)
  H_act<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist() %>% norm_vec()
  H_output<-output[4:(3+nwins)] %>% as.numeric()
  expect_equal(H_act,H_output)

  D_act<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D)
  D_act<-unlist(D_act) %>% norm_vec()
  D_output<-output[14:23]
  D_output<-as.numeric(D_output)

  expect_equal(D_act,D_output)

  h_act<-lapply(G_wins,h_stats)
  x<-rep(NA,nwins)
  haplo_act<-cbind(x,x,x,x)

  for(i in 1:nwins){
    haplo_act[i,]<-h_act[[i]]
  }

  expect_equal(as.numeric(output[24:33]),haplo_act[,1])
  expect_equal(as.numeric(output[34:43]),haplo_act[,2])
  expect_equal(as.numeric(output[44:53]),haplo_act[,3])
  expect_equal(as.numeric(output[54:63]),haplo_act[,4])
})

