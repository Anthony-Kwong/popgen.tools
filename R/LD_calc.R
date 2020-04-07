#' LD_calc function
#' 
#' Computes three LD statistics for a genome matrix. Function computes the
#' standardized LD coefficient D' for each column pair and returns the average 
#' and the maximum values of D' for the input block. D' is denoted LD in the code
#' to avoid confusion with Tajima's D. Function also computes w_max.
#' 
#' For underlying theory see
#' w_max: DOI: 10.1534/genetics.103.025387  
#' D': Lewontin RC. The Interaction of Selection and Linkage. I. General Considerations; Heterotic Models. Genetics. 1964 Jan;49(1):49-67. PMID: 17248194; PMCID: PMC1210557.
#'
#' @param G : A NumericMatrix designating a binary genome matrix consisting of 1's and 0's. 
#' Each column is a SNP. Each row is a sampled individual.  
#' 
#'
#' @return A tibble of 3 elements (LD_avg,LD_max,w_max)
#' @export
#'
#' @examples   seq <-matrix(sample(0:1, size =16 , replace = TRUE), nc = 4)
#'  LD_calc(seq)
#' @importFrom genetics LD
#' @importFrom tibble tibble
LD_calc = function(G){
  
  genotypes = matrix2genotype(G)
  #print(genotypes)
  ngeno = ncol(genotypes)
  data = genetics::LD(genotypes)
  
  # compute standardised D ----
  D = data$`D'`
  D_values = D[upper.tri(D)]
  LD_avg = mean (D_values)
  LD_max = max(D_values)
  
  #compute w_max ----
  
  r = data$`R^2`
  x = seq(2,ngeno-2)
  w = sapply(x, LD_w, r=r)
  w_max= max(w)
  
  df = tibble::tibble("LD_avg"= LD_avg, 
                      "LD_max" = LD_max,
                      "w_max" = w_max)
  #print("done")
  return(df)
}
