#' @name shrimp.snp.dat.pools
#' @docType data
#' @keywords datasets
#' @title SNP data for pools - 'pooling by phenotype' example
#' @description SNP intensity and genotype call data for DNA pools - simulated data (refer to Hamilton et al. in prep for details)
#' @usage shrimp.snp.dat.pools
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.snp.dat.pools.rda")  
#shrimp.snp.dat.pools <- sim.snp.dat.pools

shrimp.snp.dat.pools$SAMPLE_ID     <- as.integer(shrimp.snp.dat.pools$SAMPLE_ID)
shrimp.snp.dat.pools$SNP_ID     <- as.character(shrimp.snp.dat.pools$SNP_ID)
shrimp.snp.dat.pools$INTENSITY_A   <- as.numeric(shrimp.snp.dat.pools$INTENSITY_A)
shrimp.snp.dat.pools$INTENSITY_B   <- as.numeric(shrimp.snp.dat.pools$INTENSITY_B)
shrimp.snp.dat.pools$GENOTYPE     <- as.character(shrimp.snp.dat.pools$GENOTYPE)
shrimp.snp.dat.pools <- shrimp.snp.dat.pools[,c('SAMPLE_ID', 'SNP_ID', 'INTENSITY_A', 'INTENSITY_B', 'GENOTYPE')]

usethis::use_data(shrimp.snp.dat.pools, overwrite = TRUE)
