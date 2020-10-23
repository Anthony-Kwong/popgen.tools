#' ancient_sum_stats function
#' 
#' Takes in a single simulation object, ages the DNA and outputs summary statistics with the 
#' sum_stats function. The aging process consists of simulating missing information, pseduo-
#' haplotypes, deamination/DNA damage and ascertainment bias. Multiple pairs of individuals 
#' can be used for ascertainment bias. In that case, we will subset the columns of G which are
#' heterozygous in any of the pairs of individuals. See the age_DNA, pseudo_hap and
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
#' @param ascertain_indices: A list of 2-dimensional vectors. Each vector contains the indices of 2 
#' rows in G for doing the ascertainment bias. 
#' @param impute_method: A string indicating the imputation method for missingness. Options are
#' "random" and "zero." See documentation on random_impute and zero_impute for more information.
#' @param seed: Optional. A random seed for aging DNA.
#' @param denoise_method: A method for denoising the genome matrix for the purposes of computing the 
#' haplotype statistics. Default is "none". Options are "cluster" and "majority_flip". See majority_flip
#' and clus_hstats documentation for more information.
#' @param max_clus: A number between 0 and 1. Specifies the max number of clusters to consider as a fraction
#' of the total number of rows in the genome matrix. Only used with the "cluster" denoise method. Default
#' value is 0.2. 
#' 
#' @import magrittr
#' @return A one row dataframe summary stats and information for the input simulation object.
#' @export
#'
#' @examples ancient_sum_stats(sim_obj)
ancient_sum_stats <- function(sim,nwins=1,split_type="base",
                              ID,trim_sim=F,snp = NA,
                              missing_rate, trans_prop = 0.776, dmg_rate = 0.05, ascertain_indices,
                              seed = NA, impute_method, denoise_method = "none", 
                              max_clus = 0.2){
 
  #check arguments are entered correctly ----
  
  #sim
  if(class(sim)!="sim_obj"){
    stop("Argument sim is not of the class sim_obj")
  }
  
  #nwins
  if(class(nwins)!="numeric"){
    stop("Argument nwins must be a numeric")
  }
  if(floor(nwins)!=nwins){
    stop("Argument nwins must be an integer")
  }
  if(nwins<1){
    stop("Argument nwins must be a positive integer")
  }
  
  #split_type
  valid_splits=c("base","mut")
  check=split_type %in% valid_splits
  if(check!=T){
    stop("Invalid argument for split_type. Options are \"base\" or \"mut\" ")
  }
  
  #ID
  if(floor(ID)!=ID){
    stop("ID must be an integer")
  }
  if(is.numeric(ID)==F){
    stop("ID must be numeric")
  }
  
  #trim_sim
  if(trim_sim && is.na(snp)){
    stop("snp argument must an integer if trim_sim is TRUE. snp is currently NA.")
  }
  
  #snp
  if(is.na(snp)==F){
    if(floor(snp)!=snp){
      stop("ID must be an integer")
    }
    if(class(snp)!="numeric"){
      stop("ID must be numeric")
    }
  }
  
  #denoise method
  valid_impute=c("random","zero")
  check = impute_method %in% valid_impute
  if(check!=TRUE){
    stop("Invalid impute_method. Options are \"random\" or \"zero\" ")
    
  }
  
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
  
  #ascertainment bias ----
  
  #simulate ascertainment bias using all the outgroup pairs
  num_pairs <- length(ascertain_indices)
  #asc_sites_list <- vector(mode = "list", length = num_pairs)
  
  #find the indices of all the columns that are heterozygous in any of the outgroup pairs
  asc_sites_list <- lapply(seq(num_pairs), function(d){ascertain_bias(raw_G, ascertain_indices[[d]])[[2]] } )
  asc_sites <- asc_sites_list %>%
    unlist() %>% 
    unique() %>%
    sort() #sort removes any NAs which would indicate that no het sites were found in some outgroup pairs
  
  #find the indices of all the rows in the outgroup. These rows are not included for the
  #purposes of computing the summary statistics.
  asc_rows <- ascertain_indices %>%
    unlist()
  asc_G = raw_G[-c(asc_rows),asc_sites]
  
  #add missingness and deamination
  dmg_G <- age_DNA(G = asc_G, missing_rate = missing_rate,trans_prop = trans_prop,dmg_rate = dmg_rate, seed = seed)
  dmg_G <- rm_nonpoly_cols(dmg_G)[[1]]
  
  #imputation
  if(impute_method == "random"){
    set.seed(seed) #set random seed for random impute
    imp_G <- random_impute(dmg_G)
  } else if (impute_method == "zero"){
    imp_G <- zero_impute(dmg_G)
  } else {
    stop("Invalid impute_method.")
  }

  #compute the position vector for the aged genome matrix
  imp_pos <- sim$pos[asc_sites]
  
  #Before computing SS, we remove all non-polymorphic sites
  rm_G = rm_nonpoly_cols(imp_G)
  final_G = rm_G[[1]]
  final_pos = imp_pos[c(rm_G[[2]])]
  
  #Split genome matrix into windows ----
  
  #to equalise the number of SNPs across the simulations, we keep the central k SNPs around the selected mutation.
  
  #in our current pipeline, the selected mutation is always at 0.5. Later on we may change this. 
  mutation_pos = 0.5
  
  if(ncol(final_G) > snp && trim_sim){
    #find closest SNP to the selected mutation
    raw_pos<-final_pos
    snp_dist<-abs(raw_pos-mutation_pos)
    center<-which.min(snp_dist)
    
    #trim the genome matrix and pos vector
    G <- window_trim(final_G,cen=center,k=floor(snp/2))
    pos_vec<-vector_trim(final_pos, cen=center, k=floor(snp/2)) 
    
  } else {
    #no trimming of genome matrix and pos vector
    G <- final_G
    pos_vec <- final_pos
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
    
    #return NULL result if nwins > ncol(G)
    if(is.null(win_list)){
      warning("Warning: Requested number of subwindows exceed the number
              of columns in G")
      return (NA)
    }
    
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
  
  #denoise the genome matrix of pseudo haplotypes
  #if cluster if chosen, it will be done later
  if (denoise_method == "majority_flip"){
    pseudo_G = majority_flip(pseudo_G)
  }
  
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
  
  #For the cluster method, we cluster the haplotypes for computing hstats
  if(denoise_method == "cluster"){
    win_clus_vec = lapply(hap_win_list, function(M){clus_hap(M, max_clus = round( nrow(M)*max_clus) ) })
    h_values = lapply(win_clus_vec,clus_hstats)
  } else {
    h_values <- lapply(hap_win_list,h_stats)
  }
  
  
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
                             missing_rate, trans_prop, dmg_rate,
                             impute_method, denoise_method)
  
  wide_df<-cbind(sim_info,stats,base_lengths, snp_lengths)
  
  return(wide_df)

}