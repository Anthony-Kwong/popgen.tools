#' sum_stats function
#' 
#' Turns a simulation object into a data frame. The genome matrix is partitioned into subwindows and summary statistics are computed. Each subwindow is a row on the output tibble. 
#'
#' @param sim: a simulation object
#' @param nwins: number of subwindows to split each genome matrix within the simulation.
#' @param split_type: Method of splitting the genome matrix. Valid options are "base" and "mut".
#' @param ID: an ID value to group subwindows under the simulation it came from. 
#' @param snp: number of snps to include per simulation
#' @param form: output data frame in "wide" or "tall" form. This is an optional argument with default="wide". 
#' @param fun: option to apply a function over the ss across a simulation. Default is "none". The haplotype statistics (h_stats) don't get transformed. Options include "norm"
#' @importFrom tibble enframe as_tibble tibble
#' @importFrom purrr map2 pmap
#' @importFrom dplyr bind_cols
#' @import magrittr
#' @return list of summary stats for the input genome matrix
#' @export
#' @examples sum_stats(win_list)
#' This is meant to be a hidden function. Hide in final version. 
sum_stats<-function(sim,nwins,split_type,ID,snp,form="wide",fun="none"){
  
  # ensure arguments are entered correctly ----
  
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
    stop("ID must be numeric yada")
  }
  
  #snp
  if(floor(snp)!=snp){
    stop("ID must be an integer")
  }
  if(class(snp)!="numeric"){
    stop("ID must be numeric")
  }
  
  #form
  if(form!="tall" && form!="wide"){
    stop("Invalid argument:form. form must be \"wide\" or \"tall\"")
  }
  
  #fun
  if(fun!="none" && fun!="norm"){
    stop("Invalid argument:fun. See documentation for valid options")
  }
  
  
  #Quick way to see where the simulations are up to. 
  print(ID)
  
  #reject cases where there are more subwindows than SNPs. 
  if(nwins>sim$num_seg){
    txt<-paste("reject ",ID)
    print(txt)
    stop("Insufficient SNPs.",nwins, " subwindows requested and the simulation had ", sim$num_seg, " SNPs.")
  }


  
  #extract useful information from the simulation ----
  
  #if selection coefficient s=0, it is a neutral simulation
  sweep <- sim$sweep
  s_coef <- sim$s
  bottle_time1 <- sim$bottle_time1
  bottle_time2 <- sim$bottle_time2
  bottle_size1 <- sim$bottle_size1
  bottle_size2 <- sim$bottle_size2
  

  #Split genome matrix into subwindows----
  
  #to equalise the number of SNPs across the simulations, we keep the central k SNPs around the selected mutation.
  
  #in our current pipeline, the selected mutation is always at 0.5. Later on we may change this. 
  mutation_pos<-0.5
  
  if(sim$num_seg>snp){
    #find closest SNP to the selected mutation
    raw_pos<-sim$pos
    snp_dist<-abs(raw_pos-mutation_pos)
    center<-which.min(snp_dist)
    
    #trim the genome matrix and pos vector
    G<-window_trim(sim$genomes,cen=center,k=floor(snp/2))
    pos_vec<-vector_trim(raw_pos, cen=center, k=floor(snp/2)) 
    
  } else {
    #no trimming of genome matrix and pos vector
    G<-sim$genomes
    pos_vec <- sim$pos
  }
  
  if(split_type=="base"){
    win_list=winsplit_base(G,pos_vec,nwins)
  } else if (split_type=="mut"){
    win_list= sub_win(G,nwins)
  } else {
    stop("Invalid argument for split_type. Valid options are \"base\" and \"mut\".")
  }
  
  
  #Compute SS on subwindows----

  #list of basic summary statistics functions to use on the windows.These form the basis for other summary stats. 
  ss<-list(theta_h,theta_t,theta_w,var_taj)
  basic_values<-lapply(ss, function(f) sapply(win_list, function(d) f(d) ) )
  names(basic_values)<-c("theta_h","theta_t","theta_w","var_taj")
  
  #Fay and Wu H, fwh(t_t,t_h)----
  H<-purrr::map2(basic_values$theta_t,basic_values$theta_h,fwh) %>% unlist()
  #H<-H %>% tibble::enframe(name=NULL,value="FW_H") 
  
  #Tajima'D taj_D(t_t, t_w, var_taj)----
  D<-purrr::pmap(list(basic_values$theta_t,basic_values$theta_w,basic_values$var_taj),taj_D) %>% unlist()
  #D<-D %>% tibble::enframe(name=NULL,value="Taj_D") 
  
  #collect stats
  win_stats<-list("H"=H,"D"=D)
  
  #output raw values
  #win_stats<-lapply(win_stats,norm_vec)
  
  
  # Haplotype stats (h1,h2,h12,h123)----
  h_values<-lapply(win_list,h_stats) 
  names(h_values)<-string_labels("subwindow",length(h_values))
  
  #apply function over stats, note that h_stats are not transformed
  if(fun=="norm"){
    win_stats<-lapply(win_stats,norm_vec)
  }


  #tying everything together into a wide dataframe
  
  if(form=="wide"){
    
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
    final_stats<-list(win_stats,h_list) %>% unlist(recursive = FALSE)

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
    
    stats<-as.data.frame(t(df)) 
    temp<-tibble::tibble(ID,sweep,s_coef,bottle_time1,bottle_size1,
                         bottle_time2,bottle_size2)
    wide_df<-cbind(temp,stats)

    return(wide_df)
  }
  
  
  
  #tall form computation----
  
  #This tall form implementation is outdated!
  
  
  if(form=="tall"){
    #Compute distances
    
    # #Quantify distance between each subwindow and the selected mutation. Take the chromosome distance between middle of subwindow and mutation. 
    # snp_pos<-vec_split(sim$pos,nwins)
    # 
    # #preallocate memory
    # dist<-rep(NA,nwins)
    # 
    # #store distances between each subwindow and the selected mutation
    # for(i in 1:nwins){
    #   sub_win_mid<-snp_pos[i] %>% unlist() %>% median()
    #   dist[i]<-abs(sub_win_mid-mutation_pos)
    # }
    # 
    # s<-rep(s_coef,nwins)
    # sweep<-rep(sweep,nwins)
    # 
    # ID<-rep(ID,nwins)
    # 
    # #getting h_stats into the right form
    # h_df<-h_df %>% as.matrix() %>% t()
    # colnames(h_df)<-c("h1","h2","h12","h123")
    # #normalisation step
    # h_df<-apply(h_df,2,norm_vec)
    # h_df<-h_df %>% tibble::as_tibble()
    # 
    # #tying everything back together. ----
    # #Tibble is great in giving neat names without the $. 
    # # D<-norm_vec(D)
    # # H<-norm_vec(H)
    # pi_est<-basic_values$theta_t
    # df<-tibble::tibble(sweep,ID,s,dist,D,H,pi_est) %>% tibble::as_tibble()
    # tall_df<-dplyr::bind_cols(df,h_df)
    # return(tall_df)
    
    stop("tall form implementation is outdated.")
  }
  
  # change sweep and position into factors.
  # df[,c("sweep","position")]<-lapply(df[,c("sweep","position")],as.factor)

}

#inserting for testing purposes. Will remove. This bit causes trouble if left in when building. 
 # data<-readRDS("~/work/MPhil/data/hard.rds")
 # sim<-data[[70]]
 # df<-read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/toy_df.csv")
 #test<-generate_df(df,2)
 # 
 # generate_df(data,10)
 # 
 #  sim<-data[[213]]
 #  nwins=10
 #  test1<-sum_stats(sim,nwins,100)

##### This version works hurray!

# f1 <- function(x) {x}
# f2 <- function(x) {2*x}
# f3 <- function(x) {3*x}
# funs<-list(f1,f2,f3)
# 
# arg<-list(1,2,3)
# lapply(funs, function(f) sapply(arg, function(d) f(d) ) )