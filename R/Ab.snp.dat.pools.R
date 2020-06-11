#' @name ab.snp.dat.pools
#' @docType data
#' @keywords datasets
#' @title SNP data for pools - 'pooling for individual parentage assignment' example
#' @description SNP intensity and genotype call data for DNA pools - simulated data (refer to Hamilton et al. in prep for details)
#' @usage ab.snp.dat.pools
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/ab.snp.dat.pools.rda")  
#ab.snp.dat.pools <- sim.snp.dat.pools

ab.snp.dat.pools$SAMPLE_ID     <- as.integer(ab.snp.dat.pools$SAMPLE_ID)
ab.snp.dat.pools$SNP_ID     <- as.character(ab.snp.dat.pools$SNP_ID)
ab.snp.dat.pools$INTENSITY_A   <- as.numeric(ab.snp.dat.pools$INTENSITY_A)
ab.snp.dat.pools$INTENSITY_B   <- as.numeric(ab.snp.dat.pools$INTENSITY_B)
ab.snp.dat.pools$GENOTYPE     <- as.character(ab.snp.dat.pools$GENOTYPE)
ab.snp.dat.pools <- ab.snp.dat.pools[,c('SAMPLE_ID', 'SNP_ID', 'INTENSITY_A', 'INTENSITY_B', 'GENOTYPE')]

usethis::use_data(ab.snp.dat.pools, overwrite = TRUE)
