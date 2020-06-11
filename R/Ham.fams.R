#' @name Ham.fams
#' @docType data
#' @keywords datasets
#' @title Family data for small worked example in Hamilton et al. in prep
#' @description Identifies families and their sires and dams
#' @usage Ham.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

Ham.fams <- data.frame(FAMILY_ID = c(101201, 101202, 102203, 102204),
                       SIRE_ID = c(201, 202, 203, 204),
                       DAM_ID = c(101, 101, 102, 102))

usethis::use_data(Ham.fams, overwrite = TRUE)

