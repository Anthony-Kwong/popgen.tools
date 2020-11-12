#' ancient_sum_stats function
#' 
#' Takes in a single simulation object, ages the DNA and outputs the SFS for each window. 
#' See the documentation of ancient_sum_stats to see how DNA ageging is simulated. 
#'
#' @param sim: a simulation object
#' @param nwins: number of subwindows to split each genome matrix within the simulation. Default is 1.
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
#' 
#' @import magrittr
#' @return A list of NumericVectors containing the SFS for each window/block.
#' @export
#'
#' @examples ancient_sum_stats(sim_obj)
compute_sfs <- function(sim,nwins=1,missing_rate, trans_prop = 0.776, dmg_rate = 0.05, ascertain_indices,
                        seed = NA,impute_method,trim_sim=F,snp = NA){
  
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
  
  #extract useful information from the simulation ----
  
  if(is.na(seed)){
    seed = sample.int(.Machine$integer.max, size = 1)
  }
  
  #extract genome matrix
  raw_G <- sim$genomes
  
  #DNA aging component ----
  
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
  
  #add pseudo-haplodisation
  set.seed(seed)
  pseudo_G = pseudo_hap(dmg_G, seed = seed)
  
  #imputation
  if(impute_method == "random"){
    set.seed(seed) #set random seed for random impute
    imp_G <- random_impute(pseudo_G)
  } else if (impute_method == "zero"){
    imp_G <- zero_impute(pseudo_G)
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
  
  #split windows based by ~equal SNPs
  split_wins = sub_win(G,nwins)
  win_list= split_wins$windows
  
  #Compute SFS on the subwindows----
  SFS = lapply(win_list, function(M){colSums(M)})
  return (SFS)
}