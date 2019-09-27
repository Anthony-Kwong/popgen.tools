#' discoal_sim function
#'
#' @param mu mutation rate per base per generation. Warning! Check this again later in the discoal manual!
#' @param recomb_rate The recombination rate per base per generation
#' @param Ne The effective population size
#' @param genome_length The number of bases to simulate for each sample. 
#' @param samplesize Number of samples to take from the population
#' @param s selection coefficient for the selected mutation
#' @param discoal_path path to your discoal program
#' @param fix_generation number of generations ago when the selected mutation was fixed
#' @param sweep the kind of selective sweep. Input "hard" or "neutral"
#' @param seed vector of 2 numbers used for the simulations
#'
#'
#' @return an object of class sim_obj. Here are the features. cmd is the command. Seeds: the seeds used in the discoal simulation.
#' num_seg: number of segregating sites in the sampled population. pos: vector of the positions of every seg site (infinite sites model)
#' sweep: the kind of selective sweep modelled. s: the selection coefficient
#'
#' @export
#'
#' @examples sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_generation=fix,sweep=sweep)
#'
#' @importFrom stringr str_extract_all
#' @import magrittr



discoal_sim<-function(mu,recomb_rate,Ne,genome_length,samplesize,s=0,discoal_path,fix_generation,seed,sweep){
  
  #====================================================================================
  
  "Explanation of how the simulator works. 
  
  In coalescent simulations, we don't simulate individual bases but a region [0,1].
  We define length this region represents using....
  
  The nSites param within discoal technically don't represent the number of bases. They are
  points where recombination can occur. The maximum nSites for discoal is 200,000 by default. 
  
  The genome_length param in this function represents the number of base pairs to simulate. 
  To account for this, we scale up the mutation and recombination rates by genome_length. 
  
  Eg. Suppose you want to simulate a 1Mb region. The population size is N and the mutation
  rate is 1e-8 per base. 
  
  theta=4*N*(1e-8)*(1*6)
  
  The same logic applies for the scaled recombination rate rho. 
  
  When a mutation happens, we randomly select a number from U~[0,1]. This 
  will be the position of the mutation. 
  
  
  "
  
  #====================================================================================

  #setting up params for discoal command. Discoals has to scale mutation, recombination rates and selection coefficient by Ne.
  #source: Kern 2017 "discoal-a coalescent simulator with selection"

  #These params are not scaled by the number of sites.
  
  #scaled strength of selection.
  #Eg. A mutation with s=0.01 has 1% more probability of being passed down. The size of the region doesn't change this. 
  alpha=(2*Ne*s) %>% no_scientific() 
  
  #Scaled time of fixation. 
  #The time is defined by the number of coalescent units and only scales with pop size.
  tau= (fix_generation/(4*Ne)) %>% no_scientific() 
  
  #These params are scaled by the number of sites. 
  #Consider that both mutation rates and recombination rates are defined as P(event)/base.
  
  #scaled mutation rate.
  theta=(4*Ne*mu*genome_length) %>% no_scientific()  
  
  #scaled recombination rate. 
  rho=(4*Ne*recomb_rate*genome_length) %>% no_scientific() 

  #To make the code easier, we will just do one simulation at a time for now. Modify at a later date.
  nrep=1
  
  #The maximum number of  break points (nSites) of discoal.
  max_breakpoint=200000;

  #generate discoal command
  if( missing(seed) ){
    cmd=paste(discoal_path, no_scientific(samplesize), nrep,no_scientific(max_breakpoint),"-t",
              theta, "-r", rho)
  } else {

    #normally we won't input seeds. This is mainly for testing purposes.
    #break points are 200,000 max. We have scaled mutation rates and recombination rates in rho and theta to account for this.  

    cmd=paste(discoal_path, no_scientific(samplesize), nrep,no_scientific(max_breakpoint),"-t",
              theta, "-r", rho,"-d", no_scientific(seed[1]), no_scientific(seed[2]))
  }

  #continue fixing this bit
  if (sweep=="hard"){
    cmd=paste(cmd,"-a", alpha,"-ws", tau)
  }

  if(sweep=="neutral"){
    cmd=paste(cmd,"-wn", tau)
  }



  #run discoal command and save output
  #discoal has the same output format as Hudson's ms. https://snoweye.github.io/phyclust/document/msdoc.pdf.
  #The command is followed by 2 random seeds.
  sim=system(cmd,intern = T)

  #extract the random seeds

  #stringr::str_extract_all(seeds, "[0-9]+")[[1]] retrieves the numbers and returns a list of 1, the + means you have a bunch of digits following
  #as.numeric then coerces it into a numeric vector of 2 elements
  #sim[2]: the second row of the sim always gives the 2 seeds

  seeds=as.numeric(stringr::str_extract_all(sim[2], "[0-9]+")[[1]])

  #extract the number of segregating sites

  #substr(sim,1,3) gives first 3 characters of each line in sim
  #sim[substr(sim,1,3)=="seg"] gives the line with "seg"
  segsites=as.numeric(gsub(pattern="segsites:",replacement ="",sim[substr(sim,1,3)=="seg"]))

  if(segsites==0){
    return("Simulation produced no segregating sites")
  }

  #extract the positions of the segregating sites. In the infinite sites model, all positions are between 0,1

  positions=gsub(pattern="positions:",replacement ="",sim[substr(sim,1,5)=="posit"])
  #coerce the positions into numeric vector
  positions=as.numeric(stringr::str_extract_all(positions, "[0-9]\\.[0-9]+")[[1]])

  #extract genome matrix

  #the genome matrix starts right after the positions line
  start=which(substr(sim,1,3)=="pos")+1
  end=length(sim)

  #small function to take a row of genome matrix (character vector) and returns numeric vector
  string_grab<-function(s) {as.numeric(strsplit(s, split="")[[1]])}
  #use string_grab to take all the rows of the genome matrix. Transponse it as discoal displays
  #each individual as a row
  genome_matrix<-t(sapply(sim[start:end],string_grab,USE.NAMES = FALSE))
  
  if(sweep=="neutral"){
    #select coeff is 0 for neutral case
    obj<-sim_obj(cmd,seeds,segsites,positions,genome_matrix,sweep,0,fix_generation)
    return(obj)
  }
  
  #for soft and hard sweeps, we need to include the selection coefficient
  obj<-sim_obj(cmd,seeds,segsites,positions,genome_matrix,sweep,s,fix_generation)
  return(obj)
#  return(list(cmd,seeds,segsites,positions,genome_matrix))
}

#sim_obj<-function(cmd,seeds,segsites,positions,genome_matrix,sweep,select_coeff,fix_generation)
#write_rds
#pluck command, purr
