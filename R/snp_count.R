#' snp_count function
#'
#' @param sim_list: a list of simulation objects
#'
#' @return A dataframe with 2 columns (sweep_type, selection_coeff, SNP_count)
#' @importFrom dplyr bind_cols
#' @import magrittr
#' @export
#' 
#' @examples snp_count(sim_list)
snp_count<-function(sim_list){
  
  obs<-length(sim_list)
  
  #preallocate memory
  sweep_type<-list(rep(NA,obs))
  s<-list(rep(NA,obs))
  SNP<-list(rep(NA,obs))
  
  for(i in 1:obs){
    sweep_type[i]<-sim_list[[i]]$sweep
    s[i]<-sim_list[[i]]$s
    SNP[i]<-sim_list[[i]]$num_seg
  }
  
  #change predictors to be the correct form
  s<-as.numeric(s)
  SNP<-as.numeric(SNP)
  sweep_type<-as.factor(sweep_type)
  
  #tying everything together
  df<-tibble::tibble(sweep_type,s,SNP)
  
  return(df)
}

#inserted for building purposes only. Comment out before building package. 

#data<-readRDS("~/work/MPhil/data/hard.rds")
