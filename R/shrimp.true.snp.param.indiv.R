#' @name shrimp.true.snp.param.indiv
#' @docType data
#' @keywords datasets
#' @title True SNP parameter data - 'pooling by phenotype' simulation example
#' @description Assumed 'true' values for SNP parameters - used for simulations presented in Hamilton et al. in prep
#' @usage shrimp.fams
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.true.snp.param.indiv.rda")  
#shrimp.fams <- fams

shrimp.true.snp.param.indiv$SNP_ID    <- as.character(shrimp.true.snp.param.indiv$SNP_ID)   
shrimp.true.snp.param.indiv$MEAN_P_AA <- as.numeric(shrimp.true.snp.param.indiv$MEAN_P_AA)
shrimp.true.snp.param.indiv$SD_P_AA   <- as.numeric(shrimp.true.snp.param.indiv$SD_P_AA)
shrimp.true.snp.param.indiv$MEAN_P_AB <- as.numeric(shrimp.true.snp.param.indiv$MEAN_P_AB)
shrimp.true.snp.param.indiv$SD_P_AB   <- as.numeric(shrimp.true.snp.param.indiv$SD_P_AB)  
shrimp.true.snp.param.indiv$MEAN_P_BB <- as.numeric(shrimp.true.snp.param.indiv$MEAN_P_BB)
shrimp.true.snp.param.indiv$SD_P_BB   <- as.numeric(shrimp.true.snp.param.indiv$SD_P_BB) 
shrimp.true.snp.param.indiv$A_ALLELE  <- as.character(shrimp.true.snp.param.indiv$A_ALLELE) 
shrimp.true.snp.param.indiv$B_ALLELE  <- as.character(shrimp.true.snp.param.indiv$B_ALLELE)

usethis::use_data(shrimp.true.snp.param.indiv, overwrite = TRUE)

