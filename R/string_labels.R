#' string_labels function
#' 
#' Generates a string vector f_1,f_2,....,f_n
#'
#' @param f: a character string
#' @param it: number of iterations
#' @return a vector of character strings of f_1,f_2,....,f_n
#' @importFrom stringr str_c
#' @export
#' @examples string_labels("Echo",5)
string_labels<-function(f,it){
  string_vec<-stringr::str_c(f,1:it,sep="_")
  return(string_vec)
}


