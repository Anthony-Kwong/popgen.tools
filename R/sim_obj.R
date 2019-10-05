#Constructor function for the simulation object. The simulation object contains a neat
#set of useful information about a discoal simulation that can be used for analysis.

#Users aren't meant to call this function themselves.

#' sim_obj function
#' 
#' Generate a simulation object to store the data produced by discoal. 
#'
#' @param cmd the command used to generate the simulation
#' @param seeds 
#' @param segsites 
#' @param positions 
#' @param genome_matrix 
#' @param sweep 
#' @param select_coeff 
#' @param fix_generation 
#'
#' @return an object called sim_obj
#' @export
#'
#' @examples This is meant to be a hidden function.
sim_obj<-function(cmd,seeds,segsites,positions,genome_matrix,sweep,select_coeff,fix_generation){

  obj<-list(
  cmd=cmd,
  seeds=seeds,
  num_seg=segsites,
  pos=positions,
  genomes=genome_matrix,
  sweep=sweep,
  s=select_coeff
  )

  class(obj)<-"sim_obj"
  return(obj)
}

#summary/print object. Wickam R packages book.


