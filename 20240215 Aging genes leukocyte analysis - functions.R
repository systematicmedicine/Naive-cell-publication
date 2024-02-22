## 20240215 Aging genes leukocyte analysis - functions.R
# Author: Cameron
# Date: 15th February 2024
#
# Miscellaneous utility functions for sue in 20240215 aging genes leukocyte analysis.R

## Define functions
assign_stars <- function(pval){
  if (pval >= 0.05) {return("")}
  else if (pval >= 0.01) {return("*")}
  else if (pval >= 0.001) {return("**")}
  else {return("***")}
}

fix_cell_name <- function(cell.name){
  cell.name.fixed <- paste(
    toupper(substr(cell.name, 1, 1)),
    substr(cell.name, 2, nchar(cell.name)),
    sep = "")
  return(cell.name.fixed)
}