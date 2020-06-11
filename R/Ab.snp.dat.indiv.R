#' @name ab.snp.dat.indiv
#' @docType data
#' @keywords datasets
#' @title SNP data for individuals - 'pooling for individual parentage assignment' example
#' @description Individual SNP intensity and genotype call data - simulated data (refer to Hamilton et al. in prep for details)
#' @usage ab.snp.dat.indiv
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/ab.snp.dat.indiv.rda")   
#ab.snp.dat.indiv <- sim.snp.dat.indiv

ab.snp.dat.indiv$SAMPLE_ID     <- as.integer(ab.snp.dat.indiv$SAMPLE_ID)
ab.snp.dat.indiv$SNP_ID     <- as.character(ab.snp.dat.indiv$SNP_ID)
ab.snp.dat.indiv$INTENSITY_A   <- as.numeric(ab.snp.dat.indiv$INTENSITY_A)
ab.snp.dat.indiv$INTENSITY_B   <- as.numeric(ab.snp.dat.indiv$INTENSITY_B)
ab.snp.dat.indiv$A_ALLELE     <- as.character(ab.snp.dat.indiv$A_ALLELE)
ab.snp.dat.indiv$B_ALLELE     <- as.character(ab.snp.dat.indiv$B_ALLELE)
ab.snp.dat.indiv$GENOTYPE     <- as.character(ab.snp.dat.indiv$GENOTYPE)
ab.snp.dat.indiv <- ab.snp.dat.indiv[,c('SAMPLE_ID', 'SNP_ID', 'INTENSITY_A', 'INTENSITY_B', 
                                        'A_ALLELE', 'B_ALLELE', 'GENOTYPE')]

usethis::use_data(ab.snp.dat.indiv, overwrite = TRUE)

