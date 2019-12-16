#' sum_stats function
#' 
#' Turns a simulation object into a data frame. The genome matrix is partitioned into subwindows and summary statistics are computed. Each subwindow is a row on the output tibble. 
#'
#' @param sim: a simulation object
#' @param win_split: number of subwindows to split each genome matrix within the simulation.
#' @param ID: an ID value to group subwindows under the simulation it came from. 
#' @param snp: number of snps to include per simulation
#' @param form: output data frame in "wide" or "tall" form. This is an optional argument with default="wide". 
#' @importFrom tibble enframe as_tibble tibble
#' @importFrom purrr map2 pmap
#' @importFrom dplyr bind_cols
#' @import magrittr
#' @return list of summary stats for the input genome matrix
#' @export
#' @examples sum_stats(win_list)
#' This is meant to be a hidden function. Hide in final version. 
sum_stats<-function(sim,win_split,ID,snp,form="wide"){
  #Quick way to see where the simulations are up to. 
  print(ID)
  
  #reject cases where there are more subwindows than SNPs. 
  if(win_split>sim$num_seg){
    txt<-paste("reject ",ID)
    print(txt)
    return (NULL)
  }
  
  #useful information from the simulation
  
  #if selection coefficient s=0, it is a neutral simulation
  s_coef<-sim$s
  if(s_coef==0){
    sweep<-"neutral"
  } else {
    sweep<-sim$sweep
  }

  #Split genome matrix into subwindows----
  
  #to equalise the number of SNPs across the simulations, we keep the central k SNPs around the selected mutation.
  
  #in our current pipeline, the selected mutation is always at 0.5. Later on we may change this. 
  mutation_pos<-0.5
  
  if(sim$num_seg>snp){
    #find closest SNP to the selected mutation
    snp_dist<-abs(sim$pos-mutation_pos)
    center<-which.min(snp_dist)
    
    #take the k/2 snps on the left and right of the center
    lower_cutoff=round(center-(snp/2))
    upper_cutoff=round(center+(snp/2))
      
    start=max(0,lower_cutoff)
    end=min(sim$num_seg,upper_cutoff)
    
    win_list<-sub_win(sim$genomes[,start:end],win_split)
  } else {
    win_list<-sub_win(sim$genomes,win_split)
  }
  
  #Compute SS on subwindows----

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
  
  #collect stats
  win_stats<-list("H"=H,"D"=D)
  
  #output raw values
  #win_stats<-lapply(win_stats,norm_vec)
  
  
  # Haplotype stats (h1,h2,h12,h123)----
  h_values<-lapply(win_list,h_stats) 
  names(h_values)<-string_labels("subwindow",length(h_values))

  #tying everything together into a wide dataframe
  
  if(form=="wide"){
    
    #h_list stores the h_stats for each subwindow. It stores 4 vectors for the statistics h1,h2,h12,h123.
    x<-rep(NA,win_split)
    h_list<-list("h1"=x,"h2"=x,"h12"=x,"h123"=x)
    num_hstat<-length(h_values[[1]])
    
    for(i in 1:num_hstat){
      for(j in 1:win_split){
        h_list[[i]][[j]]<-h_values[[j]][[i]]
      }
    }
    
    #collect all the summary stats into a single list 
    final_stats<-list(win_stats,h_list) %>% unlist(recursive = FALSE)

    df<-lapply(final_stats,c) %>% unlist()
    
    #changing the column names. We find where each summary stat starts.  
    index=which(names(df)=="H1")
    names(df)[index:(index+win_split-1)]<-string_labels("H",win_split)

    index=which(names(df)=="D1")
    names(df)[index:(index+win_split-1)]<-string_labels("D",win_split)    
        
    index=which(names(df)=="h11")
    
    names(df)[index:(index+win_split-1)]<-string_labels("h1",win_split)
    
    index=which(names(df)=="h21")
    names(df)[index:(index+win_split-1)]<-string_labels("h2",win_split)
    
    index=which(names(df)=="h121")
    names(df)[index:(index+win_split-1)]<-string_labels("h12",win_split)
    
    index=which(names(df)=="h1231")
    names(df)[index:(index+win_split-1)]<-string_labels("h123",win_split)
    
    stats<-as.data.frame(t(df)) 
    temp<-tibble::tibble(sweep,ID,s_coef)
    wide_df<-cbind(temp,stats)

    return(wide_df)
  }
  
  
  
  #tall form computation----
  
  
  if(form=="tall"){
    #Compute distances
    
    #Quantify distance between each subwindow and the selected mutation. Take the chromosome distance between middle of subwindow and mutation. 
    snp_pos<-vec_split(sim$pos,win_split)
    
    #preallocate memory
    dist<-rep(NA,win_split)
    
    #store distances between each subwindow and the selected mutation
    for(i in 1:win_split){
      sub_win_mid<-snp_pos[i] %>% unlist() %>% median()
      dist[i]<-abs(sub_win_mid-mutation_pos)
    }
    
    s<-rep(s_coef,win_split)
    sweep<-rep(sweep,win_split)
    
    ID<-rep(ID,win_split)
    
    #getting h_stats into the right form
    h_df<-h_df %>% as.matrix() %>% t()
    colnames(h_df)<-c("h1","h2","h12","h123")
    #normalisation step
    h_df<-apply(h_df,2,norm_vec)
    h_df<-h_df %>% tibble::as_tibble()
    
    #tying everything back together. ----
    #Tibble is great in giving neat names without the $. 
    # D<-norm_vec(D)
    # H<-norm_vec(H)
    pi_est<-basic_values$theta_t
    df<-tibble::tibble(sweep,ID,s,dist,D,H,pi_est) %>% tibble::as_tibble()
    tall_df<-dplyr::bind_cols(df,h_df)
    return(tall_df)
  }
  
  #change sweep and position into factors.
  #df[,c("sweep","position")]<-lapply(df[,c("sweep","position")],as.factor)

}

#inserting for testing purposes. Will remove. This bit causes trouble if left in. 
 # data<-readRDS("~/work/MPhil/data/hard.rds")
 # sim<-data[[31]]
 #test<-generate_df(df,2)
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