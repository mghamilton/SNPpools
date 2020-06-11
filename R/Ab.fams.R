#' @name ab.fams
#' @docType data
#' @keywords datasets
#' @title Family data - 'pooling for individual parentage assignment' example
#' @description Identifies families and their sires and dams - simulated data (refer to Hamilton et al. in prep for details)
#' @usage ab.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/ab.fams.rda")  
#ab.fams <- fams

usethis::use_data(ab.fams, overwrite = TRUE)

