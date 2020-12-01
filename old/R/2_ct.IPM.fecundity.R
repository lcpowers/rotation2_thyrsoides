#' Generate second required input files
#'
#' This filters the first dataset for plants that flowered and adds the browse and seed # information
#' @param input_df Output of function 1
#' @param B_rate Browsing rate/probability
#' @param seeds_mu Average number of seeds per individual


#rm(list=ls())
ct.ipm.fecundity_fun <- function(input_df, B_rate, seeds_mu){
  
   tmp <-input_df[input_df$flow == 1,]
   
   tmp$brows.t1 <- sample(c(1, 0), size = nrow(tmp), replace = TRUE, prob = c(B_rate,1-B_rate))
   
   tmp$si.t1 <- rpois(nrow(tmp), seeds_mu)
   
   return(tmp)
  
}


