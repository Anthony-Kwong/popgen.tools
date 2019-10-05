#' calc_ss function
#' 
#' Computes a variety of popgen summary stats on a genome matrix. Requires the popgen.tools package. 
#'
#' @param windows: A list of  binary genome matrices. Columns are SNPs and rows are samples. 
#' @return list of summary stats for the input genome matrix
#' @import magrittr
#' @importFrom tibble enframe
#' @importFrom purrr pmap map2
#' @export
#'
#' @examples
calc_ss<-function(win_list){
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
 # names(h_values)<-SS_labels("sim",length(h_values))
  tidyr::gather(h_values)
  names(h_values)<-c("sim1","sim2")
  
  test<- h_values %>% tibble::as_tibble() %>% as.matrix()
 # col
  
 # %>% tidyr::gather()
  
  #names(res)<-c("theta_h","theta_t","theta_w","var_taj")
  return(res)
}

#inserting for building purposes. Will remove. 
data<-readRDS("~/work/MPhil/data/toy_set.rds")
obs<-data[[1]]
sim<-obs$genomes
n_win=2
win_list<-sub_win(sim,n_win)

