#' Generate first required input files
#'
#' This function simulates the first data set to feed into models
#' @param npops Number of population per site (single value or array. If array must have exactly one value for each site to indicate number of pops per site)
#' @param sites Names or abbreviations for sites
#' @param popsize_mu Mean population size
#' @param popsize_std Population size standard deviation
#' @param year 2004 or 2005
#' @param S_rate Survival rate
#' @param F_rate Proportion of plants in pop expected to flower
#' @param n_ros

ct.ipm_fun <- function(sites, popsize_mu, popsize_std, years, S_rate, F_rate, n_ros){
  
  data_out = NULL
  for(yy in 1:length(years)){
   
  # Loop through each site
    for(ii in 1:length(sites)){
      
      # Set population size within site
      popsize <- round(rnorm(1,popsize_mu,popsize_std))
      
      # Initialize data frame as matrix
      data <- data.frame(n = paste0(1:popsize,substr(sites[ii],1,3)),
                        site = rep(sites[ii],popsize),
                        year.t = years[yy])
      
      # Number of rosette leaves in year t
      data$nl.t <- rpois(popsize, 20)
      
      # Number of rosette leaves in year t+1
      data$nl.t1 <- data$nl.t + round(runif(popsize,-5,5))
      
      # Length of longest leaf in year t
      data$ll.t <- rpois(popsize, 15) # millimeters with mean around 15
      
      # Length of longest leaf year t+1
      data$ll.t1 <- data$ll.t + round(runif(popsize,-3,3)) # millimeters
      
      # Assign 1 for survive and 0 for mortality
      data$surv = sample(c(1, 0), size = popsize, replace = TRUE, prob = c(S_rate,1-S_rate))
      
      # TO-DO: What was the size threshold? How Should this be included?
      # Assign 1 for flowered and 0 for didn't flower
      data$flow = ifelse(data$ll.t > 15,
                         sample(c(1, 0), size = popsize, replace = TRUE, prob = c(F_rate, 1-F_rate)),
                                0)
      
      # Number of rosette leaves in year t
      data$ros <- sample(c(0, 1, 2), size = popsize, replace = TRUE, prob = c(F_rate, 1-(F_rate+0.05) ,0.05))
    
      data_out <- rbind(data_out, data)
      rm(data)
    }
    
  }
  
  return(data_out)
  
}


