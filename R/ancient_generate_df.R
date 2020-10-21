#' ancient_generate_df function
#' 
#' Takes in a list of simulation objects, simulates DNA aging and generates
#' summary statistics. See the generate_df documentation for more information. 
#'
#' @param sim_list: list of simulation objects
#' @param nwins: number of subwindows desired to split each genome matrix per simulation
#' @param split_type: Method of splitting the genome matrix. Valid options are "base" and "mut". Default "base". 
#' @param trim_sim: a logical indicating whether outer snps should be removed. Default is F.
#' @param snp : Optional. Number of snps to include per simulation. Used only when trim_sim is T.
#' @param missing_rate : Probability of elements that are randomly sampled to become NA for each row in G.
#' @param trans_prop: Proportion of columns (i.e. sites) that are chosen to be transition sites.
#' Default is 0.776.
#' @param dmg_rate: Probability of a element in a transition column changing from 0 to 1, or 1 to 0. 
#' Default is 0.05.
#' @param index: A vector of indices indicating the rows of the genome matrix which represent the outgroup.
#' @param ascertain_indices: A list of 2-dimensional vectors. Each vector contains the indices of 2 
#' rows in G for doing the ascertainment bias.
#' @param impute_method: A string indicating the imputation method for missingness. Options are
#' "random" and "zero." See documentation on random_impute and zero_impute for more information.
#' @param ID: Optional. A numeric vector if ID values to label each observation according to the simulation 
#' it came from. If argument is not used, the rows will be labelled 1,2,....
#' @param denoise_method: A method for denoising the genome matrix for the purposes of computing the 
#' haplotype statistics. Default is "none". Options are "cluster" and "majority_flip". See majority_flip
#' and clus_hstats documentation for more information.
#' 
#' @return a dataframe containting the summary statistics for the list of simulation objects
#' 
#' @importFrom purrr pmap
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples generate_df(sim_list,nwins = 10, missing_rate = 0.05, index = c(99,100))
ancient_generate_df<-function(sim_list,nwins,split_type="base",trim_sim = F,snp = NA,
                              missing_rate, trans_prop = 0.776, dmg_rate = 0.05,
                              ascertain_indices, seed = NA, impute_method, ID = NA,
                              denoise_method = "none"){
  #generate a random seed if one was not given
  if(is.na(seed)){
    seed = sample(.Machine$integer.max, 1)
  }
  num_sim<-length(sim_list)
  
  #check valid ID
  
  if(is.na(ID)){
    id <- (1:num_sim)
  } else {
    id <- ID
  }
  
  if(length(id)!=num_sim){
    stop("Vector of IDs does not match the number of simulations.")
  }
  
  set.seed(seed)
  seeds <- sample(.Machine$integer.max, num_sim)

  arg_list= list(sim_list,nwins=nwins,split_type=split_type,id,
                 trim_sim = trim_sim,snp=snp,
                 missing_rate = missing_rate, trans_prop = trans_prop,
                 dmg_rate = dmg_rate, seed = seeds,impute_method = impute_method ,
                 ascertain_indices = rep(list(ascertain_indices),num_sim),
                 denoise_method = denoise_method)
  
  df_list<-purrr::pmap(arg_list,ancient_sum_stats)
  df<-dplyr::bind_rows(df_list)
  return(df)
}