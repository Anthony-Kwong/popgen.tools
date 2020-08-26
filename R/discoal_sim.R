#' discoal_sim function
#'
#' @param mu mutation rate per base per generation. Warning! Check this again later in the discoal manual!
#' @param recomb_rate The recombination rate per base per generation
#' @param Ne The effective population size
#' @param genome_length The number of bases to simulate for each sample. 
#' @param samplesize Total number of samples to take from the population
#' @param s selection coefficient for the selected mutation. Default is 0. Note that for neutral simulations s must be 0. 
#' @param discoal_path path to your discoal program
#' @param fix_time number of generations ago when the selected mutation was fixed. 
#' @param sweep the kind of selective sweep. Options are "hard", "soft", "neutral" and "neutral_fixation". 
#' @param seed vector of 2 numbers used for the simulations
#' @param start_freq Used for soft sweeps only. The mutation spreads via drift (neutral) and becomes selected only once it has reached the starting frequency. 
#' @param popsize_changes A tibble with a size and a time (measured generations) column. The size is a multiplier for the current population size. The correponding time is the time
#' of the change in generations. Current version supports 2 changes per simulation to model bottlenecks. Discoal normally uses time in units of 4Ne but this function does the conversion. 
#' @param demes Optional.A numeric integer indicating the number of population demes to simulate. Only one deme is simulated by default.
#' @param sample_dist Optional. Only usable when simulating multiple demes. A numeric integer vector indicating the number of samples to make from each deme. 
#' The sum must be equal to samplesize.
#' @param deme_join Optional. A tibble with time, deme1, deme2 columns, indicating the time (measured in generations) to join 2 particular demes. Time is numeric. Demes are numeric
#' integers indicating the index of the deme. Note that discoal uses 0 indexing for the demes. Used to join demes if they are present. Discoal normally uses time in units of 4Ne but 
#' this function does the conversion so that time is entered in generations. 
#' @return an object of class sim_obj. Here are the features. cmd is the command. Seeds: the seeds used in the discoal simulation.
#' num_seg: number of segregating sites in the sampled population. pos: vector of the positions of every seg site (infinite sites model)
#' sweep: the kind of selective sweep modelled. s: the selection coefficient
#'
#' @export
#'
#' @examples sim<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=genome_length,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep)
#'
#' @importFrom stringr str_extract_all
#' @importFrom tibble tibble
#' @importFrom tester is_positive_integer is_numeric_vector
#' @import magrittr



discoal_sim<-function(mu,recomb_rate,Ne,genome_length,samplesize,s=0,discoal_path,
                      fix_time = NA,seed,sweep,start_freq = NA,popsize_changes = NULL,
                      demes = NA, sample_dist = NA, deme_join = NULL){
  
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
  valid_sweeps = c("hard","soft","neutral","neutral_fixation")
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
  
  #fix_time is not needed for neutral simulations
  fix_absent = is.na(fix_time)
  if(fix_absent==F && sweep=="neutral"){
    msg=paste("fix_time is not used for neutral simulations.Consider simulating under neutral_fixation to condition under a mutation getting fixed at a particular timepoint")
    stop(msg)
  }
  
  #check fit_time is present for sweep simulations
  if(fix_absent==T && sweep!="neutral"){
    msg=paste("Must specifcy fix_time for sweep simulations")
    stop(msg)
  }
  
  #change sweep type to neutral if selection coefficient is 0
  if(s==0 && (sweep!="neutral") && (sweep!="neutral_fixation")){
    msg=paste("Selection coefficient is 0. Converted sweep_type to neutral")
    warning(msg)
    sweep="neutral"
  }
  
  #check deme parameters
  if(is.na(demes)==F && is_positive_integer(demes)==F){
    msg = paste("demes must be a positive integer! It is current of class ", class(demes))
    stop(msg)
  }
  
  if(is.na(sample_dist)==F && is_numeric_vector(sample_dist)==F){
    msg = paste("sample_dist must be a numeric vector! It is current of class ", class(sample_dist))
    stop(msg)
  }
  
  if(is.na(sample_dist)==F && sum(sample_dist)!=samplesize){
    msg = paste("The sum of sample_dist must be the same as samplesize. Otherwise discoal will give segmentation fault")
    stop(msg)
  }
  
  if(is.na(demes) && is.na(sample_dist) == F){
    msg = paste("You have entered sample_dist but deme argument is missing.")
    stop(msg)
  }
  
  if(is.na(demes)==F && is.na(sample_dist)){
    msg = paste("You have entered demes but sample_dist argument is missing.")
    stop(msg)
  }
  
  if(is.na(demes)==F){
    if(length(sample_dist)!= demes){
      msg = paste("The number of elements in sample_dist must match the number of demes. 
                  sample_dist is specifying the number of samples to take from each deme.")
      stop(msg)
    }
  }

  
  #check dimensions of sample_dist and number of demes
  
  #setting up params for discoal command. Discoals has to scale mutation, recombination rates and selection coefficient by Ne.
  #source: Kern 2017 "discoal-a coalescent simulator with selection"

  #These params are not scaled by the number of sites.
  
  #scaled strength of selection.
  #Eg. A mutation with s=0.01 has 1% more probability of being passed down. The size of the region doesn't change this. 
  alpha=(2*Ne*s) %>% no_scientific() 
  
  #Scaled time of fixation. 
  #The time is defined by the number of coalescent units and only scales with pop size.
  
  if(is.na(fix_time)==F){
    tau= (fix_time/(4*Ne)) %>% no_scientific() 
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
    cmd = paste(cmd,"-a", alpha,"-ws", tau, "-f", start_freq)
  }
  
  #additional arguments too add in extra demes
  if(is.na(demes)==F){
    arg = c(demes, sample_dist)
    cmd = paste(c(cmd,"-p",arg), collapse = " ")
  }
  
  #additional argument to join demes together
  if(is.null(deme_join) == F){
    nevents = nrow(deme_join)
    deme_cmd = NULL
    scaled_times = no_scientific(deme_join$time/(4*Ne))
    for (i in 1:nevents){
      print(i)
      deme_cmd = paste(deme_cmd, "-ed", scaled_times[i], 
                       deme_join$pop1[i], deme_join$pop2[i])
    }
    cmd = paste(cmd, deme_cmd, sep = "")
  }
  
  #for population size changes current implementation works with single population only. 
  pop_index=0
  
  size_cmd = NULL
  # add popsize changes
  if(is.null(popsize_changes) == F){
    nevents = ncol(popsize_changes)

    #scale times,times in discoal are in units of 4Ne, where Ne is the popsize of the reference pop.
    #Documentation calls this No.
    scaled_times = no_scientific(popsize_changes$time/(4*Ne))

    for(i in 1:nevents){
      size_cmd = paste (size_cmd,"-en", scaled_times[i], pop_index, popsize_changes$size[i])
    }
    cmd=paste(cmd,size_cmd,sep="")
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
  
  #fill in bottleneck predictors for the constant pop case. 
  #Note we currently support 2 popsize changes. 
  if(is.null(popsize_changes)==T){
    size=c(1,1)
    time=c(0,0)
    popsize_changes=tibble::tibble(size,time)
  }
  
  
  #construct sim object
  obj<-sim_obj(cmd = cmd,seeds = seeds, segsites = segsites,positions = positions,
              genome_matrix = genome_matrix,sweep = sweep,select_coeff = s,fix_time = fix_time,
              bottle_time1 = popsize_changes$time[1], bottle_size1 = popsize_changes$size[1],
              bottle_time2 = popsize_changes$time[2], bottle_size2 = popsize_changes$size[2],
              start_freq = start_freq)
  return(obj)
  
}
