#' generate_df function
#' 
#' Takes a list of simulation objects and produces a data frame. The genome matrix of each simulation is broken into a series of non-overlapping subwindows and SS are computed. Each subwindow is a row on the final dataframe. 
#'
#' @param sim_list: list of simulation objects
#' @param nwins: number of subwindows desired to split each genome matrix per simulation
#' @param split_type: Method of splitting the genome matrix. Valid options are "base" and "mut". Default "base".
#' @param snp: number of snps to include per simulation
#' @param form: Option for dataframe to be in "tall" or "wide" form. 
#' @param fun: function to apply over summary statistics (haplotype stats excluded). See sum_stats documentation for details. 
#' @param LD_downsample: logical option to downsample genome matrix for the purposes of computing LD statistics. Default is FALSE. 
#' @param ds_prop: Used with LD_downsample. Proportion of columns to downsample without replacement from the genome matrix. 
#' @param ds_seed: Optional. Numeric seed for the LD_downsample. 
#' @importFrom purrr pmap
#' @importFrom dplyr bind_rows
#' @return dataframe containing summary statistics of all the subwindows across all the simulations. 
#' @export 
#'
#' @examples generate_df(sim_list,10)
generate_df<-function(sim_list,nwins,split_type="base",snp,form="wide",fun="none",
                      LD_downsample=F, ds_prop=NA, ds_seed=NA){
  
  #generate a random seed if one was not given
  if(is.na(ds_seed)){
    ds_seed = sample(.Machine$integer.max, 1)
  }
  
  num_sim<-length(sim_list)
  id <- (1:num_sim)
  seeds <- sample(.Machine$integer.max, num_sim)
  
  arg_list= list(sim_list,nwins=nwins,split_type=split_type,
                 id,snp=snp,form=form,fun=fun,
                 LD_downsample=LD_downsample,
                 ds_prop = ds_prop,
                 ds_seed = seeds)
  
  df_list<-purrr::pmap(arg_list,sum_stats)
  df<-dplyr::bind_rows(df_list)
  
  #change sweep and position into factors.
  df[,c("sweep","ID")]<-lapply(df[,c("sweep","ID")],as.factor)
  return(df)
}