#' Generate third required input files
#'
#' This filters the first dataset for plants that flowered and adds the browse and seed # information
#' @param nseeds number of seeds in previous year
#' @param npops Number of population per site (single value or array. If array must have exactly one value for each site to indicate number of pops per site)
#' @param sites Names or abbreviations for sites
#' 
#' 

ct.ipm.seedling_fun <- function(input_df, sites, E_rates, years){
  
  data_out = NULL
  years <- unique(input_df$year)
  
  for(yy in 1:length(years)){
    
    # Loop through each site
    for(ii in 1:length(sites)){
      
      # Filter input data for year yy and site ii
      tmp <- input_df %>% 
        filter(year == years[yy] & site == sites[ii])
      
      # Sum total number of seeds there
      n_seeds <- sum(tmp$si.t1)
      
      # Set number of seedlings established based on number of seeds from last year?
      n_sdlgs = as.numeric(rbinom(n = 1, size = n_seeds, prob = E_rates[ii]))
      
      # Initialize data frame as matrix
      data <- data.frame(n = paste0(1:n_sdlgs,substr(sites[ii],1,3),"_sdlg"),
                         site = rep(sites[ii],n_sdlgs),
                         year = years[yy])
      
      # Number of rosette leaves in year t
      data$nl.t <- rpois(n_sdlgs, 3)

      # Length of longest leaf in year t
      data$ll.t <- rpois(n_sdlgs, 15) # millimeters with mean around 15
      
      # Size in year T
      data$size.t <- round(runif(n = n_sdlgs,min = 1,max = 10))
      
      data_out <- rbind(data_out, data)
      rm(data)
    }
  }
  return(data_out)
}



