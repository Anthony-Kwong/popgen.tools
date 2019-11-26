#' sum_stats function
#' 
#' Turns a simulation object into a data frame. The genome matrix is partitioned into subwindows and summary statistics are computed. Each subwindow is a row on the output tibble. 
#'
#' @param sim: a simulation object
#' @param win_split: number of subwindows to split each genome matrix within the simulation.
#' @param ID: an ID value to group subwindows under the simulation it came from. 
#' @param snp: number of snps to include per simulation
#' @importFrom tibble enframe as_tibble tibble
#' @importFrom purrr map2 pmap
#' @importFrom dplyr bind_cols
#' @import magrittr
#' @return list of summary stats for the input genome matrix
#' @export
#' @examples sum_stats(win_list)
#' This is meant to be a hidden function. Hide in final version. 
sum_stats<-function(sim,win_split,ID,snp){
  #Quick way to see where the simulations are up to. 
  print(ID)
  
  #reject cases where there are more subwindows than SNPs. 
  if(win_split>sim$num_seg){
    txt<-paste("reject ",ID)
    print(txt)
    return (NULL)
  }
  
  # #if the number of win_split>=num_seg-1, we discard it. One column subwindows aren't useful.
  # if(win_split>=(sim$num_seg-1)){
  #   txt<-paste("reject",ID)
  #   print(txt)
  #   return(NULL)
  # }
  
  #split genome matrix into equal sized windows and store as a list.
  #to equalise the number of SNPs across the simulations, we keep the central k SNPs around the selected mutation.
  
  #in our current pipeline, the selected mutation is always at 0.5. Later on we may change this. 
  mutation_pos<-0.5
  
  if(sim$num_seg>snp){
    #find closest SNP to the selected mutation
    snp_dist<-abs(sim$pos-mutation_pos)
    center<-which.min(snp_dist)
    
    #take the k/2 snps on the left and right of the center
    win_list<-sub_win(sim$genomes[,(center-snp):(center+snp)],win_split)
  } else {
    win_list<-sub_win(sim$genomes,win_split)
  }

  #list of basic summary statistics functions to use on the windows.These form the basis for other summary stats. 
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(win_list, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  #Fay and Wu H, fwh(t_w,t_h)----
  H<-purrr::map2(basic_values$theta_w,basic_values$theta_h,fwh) %>% unlist()
  #H<-H %>% tibble::enframe(name=NULL,value="FW_H") 
  
  #Tajima'D taj_D(t_t, t_w, var_taj)----
  D<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D) %>% unlist()
  #D<-D %>% tibble::enframe(name=NULL,value="Taj_D") 
  
  # Haplotype stats (h1,h2,h12,h123)----
  h_values<-lapply(win_list,h_stats) 
  names(h_values)<-string_labels("subwindow_",length(h_values))
  h_df<- h_values %>% tibble::as_tibble()
  h_df<-h_df %>% as.matrix() %>% t()
  colnames(h_df)<-c("h1","h2","h12","h123")
  h_df<-h_df %>% tibble::as_tibble()
  
  #compute distance----
  
  #Quantify distance between each subwindow and the selected mutation. Take the chromosome distance between middle of subwindow and mutation. 
  snp_pos<-vec_split(sim$pos,win_split)
  
  #preallocate memory
  dist<-rep(NA,win_split)
  
  for(i in 1:win_split){
    sub_win_mid<-snp_pos[i] %>% unlist() %>% median()
    dist[i]<-abs(sub_win_mid-mutation_pos)
  }
  
#  dist<-dist %>% tibble::tibble()
#  names(dist)<-"dist"
  
  #Various extra details about the simulation
  snp<-sim$num_seg
  snp<-rep(snp,win_split)
  
  sweep<-sim$sweep
  sweep<-rep(sweep,win_split)
  
  ID<-rep(ID,win_split)
  
  #tying everything back together. ----
  #Tibble is great in giving neat names without the $. 
  pi_est<-basic_values$theta_t
  df<-tibble::tibble(sweep,ID,snp,dist,D,H,pi_est) %>% tibble::as_tibble()
  final_df<-dplyr::bind_cols(df,h_df)

  #change sweep and position into factors.
  #df[,c("sweep","position")]<-lapply(df[,c("sweep","position")],as.factor)

  return(final_df)
}



#inserting for building purposes. Will remove. This bit causes trouble if left in. 
 data<-readRDS("~/work/MPhil/data/hard.rds")
 sim<-data[[1]]
 test<-generate_df(df,2)
 # 
 # generate_df(data,10)
 # 
 #  sim<-data[[213]]
 #  win_split=10
 #  test1<-sum_stats(sim,win_split,100)

##### This version works hurray!

# f1 <- function(x) {x}
# f2 <- function(x) {2*x}
# f3 <- function(x) {3*x}
# funs<-list(f1,f2,f3)
# 
# arg<-list(1,2,3)
# lapply(funs, function(f) sapply(arg, function(d) f(d) ) )