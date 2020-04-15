#Constructor function for the simulation object. The simulation object contains a neat
#set of useful information about a discoal simulation that can be used for analysis.

#Users aren't meant to call this function themselves.

#' sim_obj function
#' 
#' Generate a simulation object to store the data produced by discoal. 
#'
#' @param cmd The command used to generate the simulation
#' @param seeds The two random seeds used in the discoal simulation
#' @param segsites Number of segregating sites found in sampled population
#' @param positions A vector with elements between 0,1. Designates the position of each seg site. 
#' @param genome_matrix Rows are samples. Columns are mutations/segregating sites. 
#' @param sweep Type of sweep simulated. 
#' @param select_coeff S value for the strength of the sweep. 
#' @param fix_time The time of fixation of the fixed mutation. Not needed for neutral simulations.
#' Default is 0, indicating fixation at the time of sampling. 
#' @param bottle_time1 Time (in generations) of the most recent popsize change.
#' @param bottle_size1 Multiplier for the sampled/reference population size. 
#' Refers to the most recent popsize change, looking back from the present. 
#' @param bottle_time2 Time (in generations) of the second most recent popsize change.
#' @param bottle_size2 Multiplier for the sampled/reference population size. 
#' Refers to the second most recent popsize change, looking back from the present. 
#' @return an object called sim_obj
#' @export
#'
#' @examples This is meant to be a hidden function.
sim_obj<-function(cmd,seeds,segsites,positions,genome_matrix,sweep,select_coeff,fix_time=0,
                  bottle_time1=0,bottle_size1=1,bottle_time2=0,bottle_size2=1){
  #account for 1 SNPs windows
  if(dim(genome_matrix)[1]==1){
    warning("Only 1 SNP in genome matrix.")
    genome_matrix = t(genome_matrix)
  }
  
  #check inputs
  if(is_genome_matrix(genome_matrix)==F){
    stop("Invalid genome matrix given to sim_obj constructor.")
  }
  
  obj<-list(
  cmd=cmd,
  seeds=seeds,
  num_seg=segsites,
  pos=positions,
  genomes=genome_matrix,
  sweep=sweep,
  s=select_coeff,
  fix_time=fix_time,
  bottle_time1=as.numeric(bottle_time1),
  bottle_size1=as.numeric(bottle_size1),
  bottle_time2=as.numeric(bottle_time2),
  bottle_size2=as.numeric(bottle_size2)
  )

  class(obj)<-"sim_obj"
  return(obj)
}

#summary/print object. Wickam R packages book.


