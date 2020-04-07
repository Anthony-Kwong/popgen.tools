#' LD_w function
#' 
#' Computes w for some SNP position i, in a genome matrix. 
#'
#' @param r : correlation matrix showing the correlations between each SNP 
#' position, squared. 
#' @param i :index for a column in r. i must be bigger than 1 and less than the number
#' of columns in r-1.
#'
#' @return a numeric value for the omega statistic between at SNP column i. 
#' @export
#'
#' @examples   mat <- matrix(runif(64), 8, 8)
#' mat[upper.tri(mat)==F] <- 0
#' LD_w(mat,2)
LD_w = function (r,i){
  nsam = ncol(r)
  
  #terminate if there are 3 or less SNPs. We can't compute w in this case. 
  if(ncol(r)<=3){
    warning("There are 3 or less SNPs in this window. Can't compute w values.
            Returning NA.")
    return(NA)
  }
  
  if(i<=1 || i>=(nsam-1)){
    stop("index i must be between 1 and number of columns in r-1, exclusively.")
  }
  
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
