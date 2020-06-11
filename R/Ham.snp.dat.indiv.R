#' @name Ham.snp.dat.indiv
#' @docType data
#' @keywords datasets
#' @title Individual data for small worked example in Hamilton et al. in prep
#' @description Individual SNP intensity and genotype call data
#' @usage Ham.snp.dat.indiv
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

Ham.snp.dat.indiv <- read.table(file = "extdata/Ham.snp.dat.indiv.csv",  sep = ",", header = TRUE, colClasses = "character")                      

Ham.snp.dat.indiv$B_ALLELE_FREQ <- as.numeric(Ham.snp.dat.indiv$B_ALLELE_FREQ)
Ham.snp.dat.indiv$SAMPLE_ID <- as.integer(Ham.snp.dat.indiv$SAMPLE_ID)
Ham.snp.dat.indiv$GENO_ERROR <- as.logical(Ham.snp.dat.indiv$GENO_ERROR)
Ham.snp.dat.indiv$INTENSITY_A <- as.numeric(Ham.snp.dat.indiv$INTENSITY_A)
Ham.snp.dat.indiv$INTENSITY_B <- as.numeric(Ham.snp.dat.indiv$INTENSITY_B)

usethis::use_data(Ham.snp.dat.indiv, overwrite = TRUE)

