#' sum_stats function
#' 
#' Turns a simulation object into a data frame. The genome matrix is partitioned into subwindows and summary statistics are computed. Each subwindow is a row on the output tibble. 
#'
#' @param sim: a simulation object
#' @param win_split: number of subwindows to split each genome matrix within the simulation.
#' @param ID: an ID value to group subwindows under the simulation it came from. 
#' @importFrom tibble enframe as_tibble
#' @importFrom purrr map2 pmap
#' @importFrom dplyr bind_cols
#' @import magrittr
#' @return list of summary stats for the input genome matrix
#' @export
#' @examples sum_stats(win_list)
#' This is meant to be a hidden function.
sum_stats<-function(sim,win_split,ID){
  win_list<-sub_win(sim$genomes,win_split)
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(win_list, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  #Fay and Wu H, fwh(t_w,t_h)----
  H<-purrr::map2(basic_values$theta_w,basic_values$theta_h,fwh) %>% unlist()
  #H<-H %>% tibble::enframe(name=NULL,value="FW_H") 
  
  #Tajima'D taj_D(t_t, t_w, var_taj)----
  D<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D) %>% unlist()
  #D<-D %>% tibble::enframe(name=NULL,value="Taj_D") 
  
  # Haplotype stats (h1,h2,h12)----
  h_values<-lapply(win_list,h_stats) 
  names(h_values)<-string_labels("sim",length(h_values))
  h_df<- h_values %>% tibble::as_tibble()
  h_df<-h_df %>% as.matrix() %>% t()
  colnames(h_df)<-c("h1","h2","h12")
  h_df<-h_df %>% tibble::as_tibble()
  
  #Classify position of subwindows----
  
  #The subwindow with the mutation is "mut". Immediately adjacent subwindows are "close". Otherwise "far". 
  pos<-sim$pos
  #in our current pipeline, the selected mutation is always at 0.5. Later on we may change this. 
  mutation_pos<-0.5
  sweep_pos<-which.min(abs(pos-mutation_pos))
  
  #initialise position vector
  position<-rep("far",win_split)
  
  #finding which subwindow the mutation lies in. 
  whole_win_width<-length(pos)
  width<-floor(whole_win_width/win_split)
  mut_index<-ceiling(sweep_pos/width)
  
  #updating position vector
  position[mut_index]<-"mut"
  
  if((mut_index-1)>0){
    position[mut_index-1]="close"
  }
  
  if((mut_index+1)<whole_win_width){
    position[mut_index+1]="close"
  }
  
  #position<-position %>% tibble::tibble()
  #names(position)<-"position"
  #position$position<- position$position %>% as.factor()
  
  #Various extra details about the simulation
  snp<-sim$num_seg
  snp<-rep(snp,win_split)
  
  sweep<-sim$sweep
  sweep<-rep(sweep,win_split)
  
  ID<-rep(ID,win_split)
  
  #tying everything back together. ----
  #Tibble is great in giving neat names without the $. 
  pi_est<-basic_values$theta_t
  df<-tibble::tibble(sweep,ID,snp,position,D,H,pi_est)
  df<-dplyr::bind_cols(df,h_df)
  
  #change sweep and position into factors.
  df[,c("sweep","position")]<-lapply(df[,c("sweep","position")],as.factor)

  return(df)
}

#inserting for building purposes. Will remove. This bit causes trouble if left in. 
 # data<-readRDS("~/work/MPhil/data/toy_set.rds")
 # sim<-data[[1]]
 # win_split=5

##### This version works hurray!

# f1 <- function(x) {x}
# f2 <- function(x) {2*x}
# f3 <- function(x) {3*x}
# funs<-list(f1,f2,f3)
# 
# arg<-list(1,2,3)
# lapply(funs, function(f) sapply(arg, function(d) f(d) ) )