#' @name shrimp.snp.dat.indiv
#' @docType data
#' @keywords datasets
#' @title SNP data for individuals - 'pooling by phenotype' example
#' @description Individual SNP intensity and genotype call data - simulated data (refer to Hamilton et al. in prep for details)
#' @usage shrimp.snp.dat.indiv
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

load(file = "data/shrimp.snp.dat.indiv.rda")   
#shrimp.snp.dat.indiv <- sim.snp.dat.indiv

shrimp.snp.dat.indiv$SAMPLE_ID     <- as.integer(shrimp.snp.dat.indiv$SAMPLE_ID)
shrimp.snp.dat.indiv$SNP_ID     <- as.character(shrimp.snp.dat.indiv$SNP_ID)
shrimp.snp.dat.indiv$INTENSITY_A   <- as.numeric(shrimp.snp.dat.indiv$INTENSITY_A)
shrimp.snp.dat.indiv$INTENSITY_B   <- as.numeric(shrimp.snp.dat.indiv$INTENSITY_B)
shrimp.snp.dat.indiv$A_ALLELE     <- as.character(shrimp.snp.dat.indiv$A_ALLELE)
shrimp.snp.dat.indiv$B_ALLELE     <- as.character(shrimp.snp.dat.indiv$B_ALLELE)
shrimp.snp.dat.indiv$GENOTYPE     <- as.character(shrimp.snp.dat.indiv$GENOTYPE)
shrimp.snp.dat.indiv <- shrimp.snp.dat.indiv[,c('SAMPLE_ID', 'SNP_ID', 'INTENSITY_A', 'INTENSITY_B', 
                                        'A_ALLELE', 'B_ALLELE', 'GENOTYPE')]

usethis::use_data(shrimp.snp.dat.indiv, overwrite = TRUE)

