#' @name shrimp.ped
#' @docType data
#' @keywords datasets
#' @title Pedigree data - 'pooling by phenotype' simulation example
#' @description Pedigree of individuals - used for simulations presented in Hamilton et al. in prep
#' @usage shrimp.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.ped.rda")  
#shrimp.fams <- fams

usethis::use_data(shrimp.ped, overwrite = TRUE)

