#' @name Ham.fam.set.combns
#' @docType data
#' @keywords datasets
#' @title Family set combination data for small worked example in Hamilton et al. in prep
#' @description Identifies family set combinations, family sets within combinations and families with family sets
#' @usage Ham.fam.set.combns
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

Ham.fam.set.combns <- data.frame(FAM_SET_COMBN_ID = c(1,1,1,1),
                                 FAM_SET_ID = c(1,1,2,2),
                       FAMILY_ID = c(101201, 101202, 102203, 102204))

usethis::use_data(Ham.fam.set.combns, overwrite = TRUE)

