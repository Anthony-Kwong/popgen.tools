#' discoal_sim function
#'
#' @param mu mutation rate per base per generation. Warning! Check this again later in the discoal manual!
#' @param recomb_rate The recombination rate per base per generation
#' @param Ne The effective population size
#' @param genome_length The number of bases to simulate for each sample. 
#' @param samplesize Number of samples to take from the population
#' @param s selection coefficient for the selected mutation. Default is 0. Note that for neutral simulations s must be 0. 
#' @param discoal_path path to your discoal program
#' @param fix_generation number of generations ago when the selected mutation was fixed. Default value is 0. 
#' @param sweep the kind of selective sweep. Options are "hard", "soft", "neutral" and "neutral_fixation". 
#' @param seed vector of 2 numbers used for the simulations
#' @param start_freq Used for soft sweeps only. The mutation spreads via drift (neutral) and becomes selected only once it has reached the starting frequency. 
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



discoal_sim<-function(mu,recomb_rate,Ne,genome_length,samplesize,s=0,discoal_path,fix_generation=NA,seed,sweep,start_freq=NA){
  
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

  #Check that parameters have been entered correctly
  
  #check that s is 0 for neutral simulations
  if((sweep=="neutral"||sweep=="neutral_fixation")&&s!=0){
    msg=paste("Selection coefficient must be 0 for neutral simulations. You have input s =",s, " and sweep=",sweep)
    stop(msg)
  }
  
  if(sweep=="soft" && is.na(start_freq)){
    msg=paste("start_freq must be specified for soft sweeps")
    stop(msg)
  }
  
  
  #check that a valid sweep type has been entered
  valid_sweeps=c("hard","soft","neutral","neutral_fixation")
  if((sweep %in% valid_sweeps)==F){
    msg=paste("Invalid sweep parameter. You entered sweep =",sweep," See documentation on discoal_sim.")
    stop(msg)
  }
  
  #check that starting frequency is between 0 and 1, should it be entered.
  if(is.na(start_freq)==F){
    if(start_freq<0 || start_freq>1){
      msg=paste("start_freq must be between 0 and 1")
      stop(msg)
    }
  }
  
  #fix_generation is not needed for neutral simulations
  if(is.na(fix_generation)==F && sweep=="neutral"){
    msg=paste("fix_generation is not used for neutral simulations.Consider simulating under neutral_fixation to condition under a mutation getting fixed at a particular timepoint")
    stop(msg)
  }
  
  #change sweep type to neutral if selection coefficient is 0
  if(s==0 && (sweep!="neutral") && (sweep!="neutral_fixation")){
    msg=paste("Selection coefficient is 0. Converted sweep_type to neutral")
    warning(msg)
    sweep="neutral"
  }
  
  #warning when s is unspecified for selective sweeps
  
  
  #setting up params for discoal command. Discoals has to scale mutation, recombination rates and selection coefficient by Ne.
  #source: Kern 2017 "discoal-a coalescent simulator with selection"

  #These params are not scaled by the number of sites.
  
  #scaled strength of selection.
  #Eg. A mutation with s=0.01 has 1% more probability of being passed down. The size of the region doesn't change this. 
  alpha=(2*Ne*s) %>% no_scientific() 
  
  #Scaled time of fixation. 
  #The time is defined by the number of coalescent units and only scales with pop size.
  
  if(is.na(fix_generation)==F){
    tau= (fix_generation/(4*Ne)) %>% no_scientific() 
  }

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
  
  #additional arguments for hard, soft and neutral_fixation
  #ordinary neutral coalescent simulations don't require further arguments. 

  if (sweep=="hard"){
    cmd=paste(cmd,"-a", alpha,"-ws", tau)
  }

  if(sweep=="neutral_fixation"){
    cmd=paste(cmd,"-wn", tau)
  }
  
  if(sweep=="soft"){
    cmd=paste(cmd,"-a", alpha,"-ws", tau, "-f", start_freq)
  }

  #print(cmd)

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
  
  if(sweep=="neutral" || sweep=="neutral_fixation"){
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
