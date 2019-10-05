#' sum_stats function
#' 
#' Computes a variety of popgen summary stats on a genome matrix list.
#'
#' @param windows: A list of  binary genome matrices. Columns are SNPs and rows are samples. 
#' @importFrom tibble enframe as_tibble
#' @importFrom purrr map2 pmap
#' @importFrom dplyr bind_cols
#' @return list of summary stats for the input genome matrix
#' @export
#' @examples sum_stats(win_list)
#' This is meant to be a hidden function.
sum_stats<-function(win_list){
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(win_list, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  #Fay and Wu H, fwh(t_w,t_h)
  H<-purrr::map2(basic_values$theta_w,basic_values$theta_h,fwh) %>% unlist()
  H<-H %>% tibble::enframe(name=NULL,value="FW_H") 
  
  #Tajima'D taj_D(t_t, t_w, var_taj)
  D<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D) %>% unlist()
  D<-D %>% tibble::enframe(name=NULL,value="Taj_D") 
  
  # Haplotype stats (h1,h2,h12)
  h_values<-lapply(win_list,h_stats) 
  names(h_values)<-string_labels("sim",length(h_values))
  h_df<- h_values %>% tibble::as_tibble()
  h_df<-h_df %>% as.matrix() %>% t()
  colnames(h_df)<-c("h1","h2","h12")
  h_df<-h_df %>% tibble::as_tibble()
  
  #tying everything back together. Tibble is great in giving neat names without the $. 
  pi_est<-tibble::tibble(basic_values$theta_t)
  names(pi_est)<-"pi"
  df<-dplyr::bind_cols(D,H,h_df,pi_est)

  return(df)
}

#inserting for building purposes. Will remove. This bit causes trouble if left in. 
# data<-readRDS("~/work/MPhil/data/toy_set.rds")
# obs<-data[[1]]
# sim<-obs$genomes
# n_win=2
# win_list<-sub_win(sim,n_win)

##### This version works hurray!

# f1 <- function(x) {x}
# f2 <- function(x) {2*x}
# f3 <- function(x) {3*x}
# funs<-list(f1,f2,f3)
# 
# arg<-list(1,2,3)
# lapply(funs, function(f) sapply(arg, function(d) f(d) ) )