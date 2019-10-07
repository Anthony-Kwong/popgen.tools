#' generate_df function
#' 
#' Takes a list of simulation objects and produces a data frame. The genome matrix of each simulation is broken into a series of non-overlapping subwindows and SS are computed. Each subwindow is a row on the final dataframe. 
#'
#' @param sim_list list of simulation objects
#' @param win_split number of subwindows desired to split each genome matrix per simulation
#' @importFrom purrr pmap
#' @importFrom dplyr bind_rows
#' @return dataframe containing summary statistics of all the subwindows across all the simulations. 
#' @export 
#'
#' @examples generate_df(sim_list,10)
generate_df<-function(sim_list,win_split){
  num_sim<-length(sim_list)
  
  #generate IDs for each simulation
  x<-(1:num_sim)
  split<-rep(win_split,num_sim)
  df_list<-purrr::pmap(list(sim_list,split,x),sum_stats)
  df<-dplyr::bind_rows(df_list)
  
  #change sweep and position into factors.
  df[,c("sweep","ID","position")]<-lapply(df[,c("sweep","ID","position")],as.factor)
  return(df)
}