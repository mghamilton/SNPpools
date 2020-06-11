#' @name ab.fam.set.combns.by.pool
#' @docType data
#' @keywords datasets
#' @title Family set combination by pool data - 'pooling for individual parentage assignment' example
#' @description Identifies family set combinations with pooled samples - simulated data (refer to Hamilton et al. in prep for details)
#' @usage ab.fam.set.combns.by.pool
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/ab.fam.set.combns.by.pool.rda")  
#ab.fam.set.combns.by.pool <- fam.set.combns.by.pool

usethis::use_data(ab.fam.set.combns.by.pool, overwrite = TRUE)

