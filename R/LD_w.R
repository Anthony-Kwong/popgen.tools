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
#' mat[upper.tri(mat)==F] <- NA
#' LD_w(mat,2)
#' @importFrom tester is_numeric_matrix is_integer has_NA
LD_w = function (r,i){
  
  #check inputs
  
  if(tester::is_numeric_matrix(r)==F){
    input_class = class(r)
    msg = paste0("Input G must be a numeric matrix.
                 Currently, G is a ", input_class)
    stop(msg)
  }
  
  #ensure the lower triangular and diagonal elements are NA
  lower_elements = c(diag(r), r[lower.tri(r)])
  check = is.na(lower_elements)
  if(all(check)==F){
    stop("Input r must have NA values on the diagonal and the lower triangular elements.
         It is a correlation matrix.")
  }
  
  nsam = ncol(r)
  
  #terminate if there are 3 or less SNPs. We can't compute w in this case. 
  if(ncol(r)<=3){
    warning("There are 3 or less SNPs in this window. Can't compute w values.
            Returning NA.")
    return(NA)
  }
  
  if(tester::is_integer(i)==F){
    stop("i must be an integer. It is a column index.")
  }
  
  if(i<=1 || i>=(nsam-1)){
    stop("index i must be between 1 and number of columns in r-1, exclusively.")
  }
  
  L = r[,1:i]
  R = r[,(i+1):nsam]
  
  #compute nominator term
  top_sums = sum(L, na.rm = T) + sum(R[(i+1):nsam,],na.rm = T)
  top_norm = choose(i,2) + choose(nsam-i,2)
  
  nominator = top_sums/top_norm
  
  #compute denominator term
  bot_sum = sum(R[(1:i),], na.rm = T)
  bot_norm = i*(nsam - i)
  
  denominator = bot_sum/bot_norm
  
  return (nominator/denominator)
}
