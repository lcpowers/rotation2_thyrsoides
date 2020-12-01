# 4. "IPM.establishment.data.txt"
#	columns:	site == Population abbreviation (here: FU == Furka Pass, SP == Schynige Platte)
#	 		plot == Number of plot
#			year.t == here: 2004 or 2005
#			seeds.t == seed production in each plot
#			sdl.t1 == seedlings in year t+1 in each plot

#' Generate fourth required input files
#'
#' This filters the first dataset for plants that flowered and adds the browse and seed # information
#' @param nseeds number of seeds in previous year
#' @param npops Number of population per site (single value or array. If array must have exactly one value for each site to indicate number of pops per site)
#' @param sites Names or abbreviations for sites
#' 
#'