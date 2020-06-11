#' @name shrimp.fams
#' @docType data
#' @keywords datasets
#' @title Family data - 'pooling by phenotype' example
#' @description Identifies families and their sires and dams - simulated data (refer to Hamilton et al. in prep for details)
#' @usage shrimp.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.fams.rda")  
#shrimp.fams <- fams

usethis::use_data(shrimp.fams, overwrite = TRUE)

