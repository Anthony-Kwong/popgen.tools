#' ancient_sum_stats function
#' 
#' Takes in a single simulation object, ages the DNA and outputs summary statistics with the 
#' sum_stats function. The aging process consists of simulating missing information, pseduo-
#' haplotypes, deamination/DNA damage and ascertainment bias. See the age_DNA, pseudo_hap and
#' ascertain_bias function documentation for more information. See sum_stats function 
#' documentation to see how summary statistics are computed. 
#'
#' @param sim: a simulation object
#' @param nwins: number of subwindows to split each genome matrix within the simulation. Default is 1.
#' @param split_type: Method of splitting the genome matrix. Valid options are "base" and "mut". Default is "base".
#' @param ID: a numeric ID value to group subwindows under the simulation it came from. 
#' @param trim_sim: a logical indicating whether outer snps should be removed. Default is F.
#' @param snp: number of snps to include per simulation
#' @param missing_rate : Probability of elements that are randomly sampled to become NA for each row in G.
#' @param trans_prop: Proportion of columns (i.e. sites) that are chosen to be transition sites.
#' Default is 0.776.
#' @param dmg_rate: Probability of a element in a transition column changing from 0 to 1, or 1 to 0. 
#' Default is 0.05.
#' @param index: A vector of indices indicating the rows of the genome matrix which represent the outgroup.
#' @param seed: Optional. A random seed for aging DNA.
#' 
#' @return A one row dataframe summary stats and information for the input simulation object.
#' @export
#'
#' @examples ancient_sum_stats(sim_obj)
ancient_sum_stats <- function(sim,nwins=1,split_type="base",
                              ID,trim_sim=F,snp = NA,
                              missing_rate, trans_prop= 0.776, dmg_rate = 0.05, index,
                              seed = NA){
 
  #extract useful information from the simulation ----
  
  sweep <- sim$sweep
  s_coef <- sim$s
  bottle_time1 <- sim$bottle_time1
  bottle_time2 <- sim$bottle_time2
  bottle_size1 <- sim$bottle_size1
  bottle_size2 <- sim$bottle_size2
  start_freq <- sim$start_freq
  
  
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1)
  }
  
  #DNA aging component ----
  raw_G <- sim$genomes
  asc_bias <- ascertain_bias(raw_G, index)
  dmg_G <- age_DNA(G = asc_bias[[1]], missing_rate = missing_rate,trans_prop = trans_prop,dmg_rate = dmg_rate, seed = seed)
  set.seed(seed)
  imp_G <- random_impute(dmg_G)
  imp_pos <- sim$pos[ asc_bias[[2]] ]
  
  #Split genome matrix into windows ----
  
  #to equalise the number of SNPs across the simulations, we keep the central k SNPs around the selected mutation.
  
  #in our current pipeline, the selected mutation is always at 0.5. Later on we may change this. 
  mutation_pos<-0.5
  
  if(ncol(imp_G) > snp && trim_sim){
    #find closest SNP to the selected mutation
    raw_pos<-imp_pos
    snp_dist<-abs(raw_pos-mutation_pos)
    center<-which.min(snp_dist)
    
    #trim the genome matrix and pos vector
    G<-window_trim(imp_G,cen=center,k=floor(snp/2))
    pos_vec<-vector_trim(imp_pos, cen=center, k=floor(snp/2)) 
    
  } else {
    #no trimming of genome matrix and pos vector
    G<-imp_G
    pos_vec <- imp_pos
  }
  
  if(split_type=="base"){
    split_wins = winsplit_base(G,pos_vec,nwins)
    win_list = split_wins$windows
    base_length = split_wins$base_length
    
    #compute number of SNPs in each block. This is more for checking purposes.
    snp_lengths = sapply(win_list, ncol)
    snp_lengths = as.data.frame(snp_lengths) %>% t()
    colnames(snp_lengths) = string_labels("block_snp_length", nwins)
    row.names(snp_lengths) = NULL
    
  } else if (split_type=="mut"){
    #split windows based by ~equal SNPs
    split_wins = sub_win(G,nwins)
    win_list= split_wins$windows
    
    #compute the length of each block in bases
    j = seq(2, nwins + 1)
    base_lengths = sapply(j, function(j){ pos_vec[ split_wins$bounds[j] ] - 
        pos_vec[ split_wins$bounds[j-1] ]})
    base_lengths = as.data.frame(base_lengths) %>% t()
    colnames(base_lengths) = string_labels("block_base_length", nwins)
    row.names(base_lengths) = NULL
    
    #compute number of SNPs in each block. This is more for checking purposes.
    snp_lengths = sapply(win_list, ncol)
    snp_lengths = as.data.frame(snp_lengths) %>% t()
    colnames(snp_lengths) = string_labels("block_snp_length", nwins)
    row.names(snp_lengths) = NULL
    
  } else {
    stop("Invalid argument for split_type. Valid options are \"base\" and \"mut\".")
  }
  
  #Compute SS on subwindows----
  
  #list of basic summary statistics functions to use on the windows.These form the basis for other summary stats. 
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(win_list, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  #Fay and Wu H, fwh(t_t,t_h) ----
  H<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  #Tajima'D taj_D(t_t, t_w, var_taj)----
  D<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D) %>% unlist()
  
  #compute haplotype stats----
  set.seed(seed)
  pseudo_G = pseudo_hap(G, seed = seed)
  
  #split the genome matrix of pseudo haplotypes
  if(split_type=="base"){
    split_wins = winsplit_base(pseudo_G,pos_vec,nwins)
    hap_win_list = split_wins$windows
    
  } else if (split_type=="mut"){
    #split windows based by ~equal SNPs
    split_wins = sub_win(pseudo_G,nwins)
    hap_win_list= split_wins$windows
    
  } else {
    stop("Invalid argument for split_type. Valid options are \"base\" and \"mut\".")
  }
  
  #skip genome matrix check because pseudo_hap can produce columns with just 0's. This is ok.
  h_values <- lapply(hap_win_list,h_stats) 
  names(h_values) <- string_labels("subwindow",length(h_values))
  
  #h_list stores the h_stats for each subwindow. It stores 4 vectors for the statistics h1,h2,h12,h123.
  x<-rep(NA,nwins)
  h_list<-list("h1"=x,"h2"=x,"h12"=x,"h123"=x)
  num_hstat<-length(h_values[[1]])
  
  for(i in 1:num_hstat){
    for(j in 1:nwins){
      h_list[[i]][[j]]<-h_values[[j]][[i]]
    }
  }
  
  #collect all the summary stats into a single list 
  final_stats<-list("D"= D,"H" = H,h_list) %>% unlist(recursive = FALSE)
  
  df<-lapply(final_stats,c) %>% unlist()
  
  #changing the column names. We find where each summary stat starts.  
  index=which(names(df)=="H1")
  names(df)[index:(index+nwins-1)]<-string_labels("H",nwins)
  
  index=which(names(df)=="D1")
  names(df)[index:(index+nwins-1)]<-string_labels("D",nwins)  
  
  index=which(names(df)=="h11")
  names(df)[index:(index+nwins-1)]<-string_labels("h1",nwins)
  
  index=which(names(df)=="h21")
  names(df)[index:(index+nwins-1)]<-string_labels("h2",nwins)
  
  index=which(names(df)=="h121")
  names(df)[index:(index+nwins-1)]<-string_labels("h12",nwins)
  
  index=which(names(df)=="h1231")
  names(df)[index:(index+nwins-1)]<-string_labels("h123",nwins)
  
  stats <- as.data.frame(t(df)) 
  sim_info <- tibble::tibble(ID,sweep,s_coef,bottle_time1,bottle_size1,
                             bottle_time2,bottle_size2,start_freq,
                             missing_rate, trans_prop, dmg_rate)
  wide_df<-cbind(sim_info,stats)
  
  if(split_type == "mut"){
    wide_df = cbind(wide_df, base_lengths, snp_lengths)
  }
  
  if(split_type == "base") {
    wide_df = cbind(wide_df, base_length, snp_lengths)
  }
  
  return(wide_df)

}