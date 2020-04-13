#' LD_calc function
#' 
#' Computes three LD statistics for a genome matrix. Function computes the
#' standardized LD coefficient D' for each column pair and returns the average 
#' and the maximum values of D' for the input block. D' is denoted LD in the code
#' to avoid confusion with Tajima's D. Function also computes w_max.
#' 
#' @references w_max: DOI: 10.1534/genetics.103.025387  
#' D': Lewontin RC. The Interaction of Selection and Linkage. I. General Considerations; Heterotic Models. Genetics. 1964 Jan;49(1):49-67. PMID: 17248194; PMCID: PMC1210557.
#' Zns: Kelly, J.K., 1997. A test of neutrality based on interlocus associations. Genetics, 146(3), pp.1197-1206. Vancouver	
#'
#' @param G : A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. 
#' Each column is a SNP. Each row is a sampled individual.  
#' 
#'
#' @return A tibble of 4 elements (LD_avg,LD_max,w_max, Zns)
#' @export
#'
#' @examples   seq <-matrix(sample(0:1, size =16 , replace = TRUE), nc = 4)
#'  LD_calc(seq)
#' @importFrom genetics LD
#' @importFrom tibble tibble
#' @importFrom tester has_NA
LD_calc = function(G){
  
  #check input ----
  if(tester::has_NA(G)){
    warning("Input matrix is NA. Returning NA for LD stats")
    df = tibble::tibble("LD_avg"= NA, 
                   "LD_max" = NA,
                   "w_max" = NA,
                   "Zns" = NA)
    return(df)
  }
  
  if(is_genome_matrix(G)==F){
    input_class = class(G)
    msg = paste0("Input G must be a valid genome matrix.
                 Currently, G is a ", input_class)
    stop(msg)
  }
  
  #compute LD stats ----
  genotypes = matrix2genotype(G)
  ngeno = ncol(genotypes)
  data = genetics::LD(genotypes)
  
  # compute standardised D ----
  # we take the absolute value as we are interested in the magnitude not the size. 
  
  D = abs(data$`D'`)
  D_values = D[upper.tri(D)]
  LD_avg = mean (D_values)
  LD_max = max(D_values)
  
  #compute w_max ----
  
  r = data$`R^2`
  x = seq(2,ngeno-2)
  w = sapply(x, LD_w, r=r)
  w_max= max(w)
  
  #compute Zns ----
  r_values = r[upper.tri(r)]
  Zns = mean(r_values)
  
  df = tibble::tibble("LD_avg"= LD_avg, 
                      "LD_max" = LD_max,
                      "w_max" = w_max,
                      "Zns" = Zns)
  #print("done")
  return(df)
}
