#' @name ab.fam.set.combns
#' @docType data
#' @keywords datasets
#' @title Family set combination data - 'pooling for individual parentage assignment' example
#' @description Identifies family set combinations, family sets within combinations and families with family sets - simulated data (refer to Hamilton et al. in prep for details)
#' @usage ab.fam.set.combns
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/ab.fam.set.combns.rda")  
#ab.fam.set.combns <- fam.set.combns

usethis::use_data(ab.fam.set.combns, overwrite = TRUE)