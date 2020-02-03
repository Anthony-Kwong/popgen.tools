#' generate_df function
#' 
#' Takes a list of simulation objects and produces a data frame. The genome matrix of each simulation is broken into a series of non-overlapping subwindows and SS are computed. Each subwindow is a row on the final dataframe. 
#'
#' @param sim_list list of simulation objects
#' @param win_split number of subwindows desired to split each genome matrix per simulation
#' @param snp: number of snps to include per simulation
#' @param form: Option for dataframe to be in "tall" or "wide" form. 
#' @param fun: function to apply over summary statistics (haplotype stats excluded). See sum_stats documentation for details. 
#' @importFrom purrr pmap
#' @importFrom dplyr bind_rows
#' @return dataframe containing summary statistics of all the subwindows across all the simulations. 
#' @export 
#'
#' @examples generate_df(sim_list,10)
generate_df<-function(sim_list,win_split,snp,form="wide",fun="none"){
  num_sim<-length(sim_list)
  
  #extend arguments into vectors for pmap
  #generate IDs for each simulation
  id<-(1:num_sim)
  split<-rep(win_split,num_sim)
  snp<-rep(snp,num_sim)
  form<-rep(form,num_sim)
  fun<-rep(fun,num_sim)
  
  df_list<-purrr::pmap(list(sim_list,split,id,snp,form,fun),sum_stats)
  df<-dplyr::bind_rows(df_list)
  
  #change sweep and position into factors.
  df[,c("sweep","ID")]<-lapply(df[,c("sweep","ID")],as.factor)
  return(df)
}