#' @name shrimp.map
#' @docType data
#' @keywords datasets
#' @title Map data - 'pooling by phenotype' simulation example
#' @description Genetic map file for each SNP - used for simulations presented in Hamilton et al. in prep
#' @usage shrimp.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.map.rda")  
#shrimp.fams <- fams

usethis::use_data(shrimp.map, overwrite = TRUE)

