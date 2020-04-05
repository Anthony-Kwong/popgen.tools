#' LD_w function
#' 
#' Computes w for some SNP position i, in a genome matrix. 
#'
#' @param r : correlation matrix showing the correlations between each SNP 
#' position, squared. 
#' @param i index for a column in r.
#'
#' @return a numeric value for the omega statistic between at SNP column i. 
#' @export
#'
#' @examples   mat <- matrix(runif(64), 8, 8)
#' mat[upper.tri(mat)==F] <- 0
#' LD_w(mat,2)
LD_w = function (r,i){
  
  nsam = ncol(r)
  r[is.na(r)] = 0
  L = r[,1:i]
  R = r[,(i+1):nsam]
  
  #compute nominator term
  top_sums = sum(L) + sum(R[(i+1):nsam,])
  top_norm = choose(i,2) + choose(nsam-i,2)
  
  nominator = top_sums/top_norm
  
  #compute denominator term
  bot_sum = sum(R[(1:i),])
  bot_norm = i*(nsam - i)
  
  denominator = bot_sum/bot_norm
  
  return (nominator/denominator)
}
