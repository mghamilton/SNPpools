#' @name sim.fam.sets
#' @docType data
#' @keywords datasets
#' @title Assignment of families to family sets - 'pooling by phenotype' simulation example
#' @description Assigns families to family sets - used for simulations presented in Hamilton et al. in prep
#' @usage shrimp.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.sim.fam.sets.rda")  
#shrimp.fams <- fams

usethis::use_data(shrimp.sim.fam.sets, overwrite = TRUE)

