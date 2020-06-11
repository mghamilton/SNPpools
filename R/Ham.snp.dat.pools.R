#' @name Ham.snp.dat.pools
#' @docType data
#' @keywords datasets
#' @title Data for DNA pools from small worked example in Hamilton et al. in prep
#' @description SNP intensity and genotype call data for DNA pools
#' @usage Ham.snp.dat.pools
#' @format Data frame
#' @author Matthew Hamilton <matthew.hamilton@csiro.au>
#' @source Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data
#' @references Hamilton et al. (in prep) Methods of parentage assignment for pooled samples using low-density SNP data

Ham.snp.dat.pools <- read.table(file = "extdata/Ham.snp.dat.pools.csv",  sep = ",", header = TRUE, colClasses = "character")                      

Ham.snp.dat.pools$B_ALLELE_FREQ <- as.numeric(Ham.snp.dat.pools$B_ALLELE_FREQ)
Ham.snp.dat.pools$SAMPLE_ID <- as.integer(Ham.snp.dat.pools$SAMPLE_ID)
Ham.snp.dat.pools$GENO_ERROR <- as.logical(Ham.snp.dat.pools$GENO_ERROR)
Ham.snp.dat.pools$INTENSITY_A <- as.numeric(Ham.snp.dat.pools$INTENSITY_A)
Ham.snp.dat.pools$INTENSITY_B <- as.numeric(Ham.snp.dat.pools$INTENSITY_B)

usethis::use_data(Ham.snp.dat.pools, overwrite = TRUE)