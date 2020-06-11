# June 2020
# Matthew Hamilton

#' snp.param.pools.fun
#' 
#' This function generates SNP parameters (specifically allelic proportion means and standard deviations) 
#' for pooled genotypes as described in Hamilton et al. (in prep).  
#' 
#' @param snp.param.indiv is the output of snp.param.indiv.fun.  That is, it is a data frame with the following headings (class in parentheses):
#' \itemize{
#'  \item{'SNP_ID' is the SNP identifier (character).} 
#'  \item{'N_AA' is the count of homozygous A (AA) genotypes (integer).}
#'  \item{'MEAN_P_AA' is the mean of allelic proportion for homozygous A genotypes (numeric).}
#'  \item{'SD_P_AA' is the standard deviation of allelic proportion for homozygous A genotypes (numeric).} 
#'  \item{'N_AB' is the count of heterozygous (AB) genotypes (integer).}
#'  \item{'MEAN_P_AB' is the mean of allelic proportion for heterozygous (AB) genotypes (numeric).}
#'  \item{'SD_P_AB' is the standard deviation of allelic proportion for heterozygous (AB) genotypes (numeric).}  
#'  \item{'N_BB' is the count of homozygous B (BB) genotypes (integer).}
#'  \item{'MEAN_P_BB' is the mean of allelic proportion for homozygous B genotypes (numeric).}
#'  \item{'SD_P_BB' is the standard deviation of allelic proportion for homozygous B genotypes (numeric).} 
#'  \item{'WELCH_A' is the welsh statistic for the interval between MEAN_P_AA and MEAN_P_AB (numeric).}
#'  \item{'WELCH_B' is the welsh statistic for the interval between MEAN_P_AB and MEAN_P_BB (numeric).}
#'  \item{'A_ALLELE_FREQ' is the A allele frequency computed from genotype counts (numeric).}
#'  \item{'B_ALLELE_FREQ' is the B allele frequency computed from genotype counts (numeric).}
#'  \item{'A_ALLELE' is the base represented by allele A (i.e. 'A', 'C', 'G' or 'T') (character).}
#'  \item{'B_ALLELE' is the base represented by allele B (i.e. 'A', 'C', 'G' or 'T') (character).}
#' } 

#' @param n.in.pools is an integer representing the number of individuals in pools.    

#' @return 1. A data frame containing estimates of SNP sepecific parameters for pooled genotypes (refer to Hamilton 2020).  
#' The number of columns in the data frame is dependent on the number of possible pooled genotypes which is determined from n.in.pools. 
#' The following is an example of output where n.in.pools = 2. 
#' \itemize{
#'  \item{'SNP_ID' is the SNP identifier.} 
#'  \item{'MEAN_P_AAAA' is the mean of allelic proportion for homozygous A genotypes.}
#'  \item{'SD_P_AAAA' is the standard deviation of allelic proportion for homozygous A genotypes.} 
#'  \item{'MEAN_P_AAAB' is the mean of allelic proportion for unordered AAAB genotypes.}
#'  \item{'SD_P_AAAB' is the standard deviation of allelic proportion for unordered AAAB genotypes.} 
#'  \item{'MEAN_P_AABB' is the mean of allelic proportion for unordered AABB genotypes.}
#'  \item{'SD_P_AABB' is the standard deviation of allelic proportion for unordered AABB genotypes.}  
#'  \item{'MEAN_P_ABBB' is the mean of allelic proportion for unordered ABBB genotypes.}
#'  \item{'SD_P_ABBB' is the standard deviation of allelic proportion for unordered ABBB genotypes.}
#'  \item{'MEAN_P_BBBB' is the mean of allelic proportion for homozygous B genotypes.}
#'  \item{'SD_P_BBBB' is the standard deviation of allelic proportion for homozygous B genotypes.}   
#'  \item{'A_ALLELE' is the base represented by allele A (i.e. 'A', 'C', 'G' or 'T').}
#'  \item{'B_ALLELE' is the base represented by allele B (i.e. 'A', 'C', 'G' or 'T').}
#' } 

#' @examples
#' #Retrieve data for small worked example from Hamilton 2020
#' data(Ham.snp.dat.indiv)
#' 
#' #Compute SNP parameters
#' Ham.snp.param.indiv <- snp.param.indiv.fun(Ham.snp.dat.indiv)
#' snp.param.pools.fun(Ham.snp.param.indiv, n.in.pools = 2)

#' @references Hamilton MG (2020) Maximum likelihood parentage assignment using quantitative genotypes

#' @export

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

snp.param.pools.fun <- function(snp.param.indiv , #from snp.gen.param.fun
                                n.in.pools
) {

  # load required packages
  if("dplyr" %in% installed.packages()[, "Package"] == FALSE) {install.packages("dplyr", repos='https://cran.csiro.au/')} 
  library(dplyr)
  
  if(sum(c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB", "A_ALLELE", "B_ALLELE") %in% 
         colnames(snp.param.indiv)) != 9) {
    stop("snp.param.indiv input must be a data frame containing the following headings: SNP_ID, MEAN_P_AA, SD_P_AA, MEAN_P_AB, SD_P_AB, MEAN_P_BB, SD_P_BB, A_ALLELE, B_ALLELE")
  }
  
  #generate alleles data frame
  alleles.by.snp <- snp.param.indiv[,c("SNP_ID","A_ALLELE", "B_ALLELE")]
  
  #Generate list of genotypes
  genotypes <- genotypes.fun(n=(2*n.in.pools))
  
  #Generate means
  means <- matrix(NA, nrow = nrow(snp.param.indiv), ncol = length(genotypes) + 1)
  means <- as.data.frame(means)
  colnames(means) <- c("SNP_ID", paste("MEAN_P_",genotypes, sep=""))
  means[,"SNP_ID"] <- snp.param.indiv[,"SNP_ID"]
  
  means[,2]                   <- snp.param.indiv[,"MEAN_P_AA"] #homozygous A
  means[,(ncol(means)/2 + 1)] <- snp.param.indiv[,"MEAN_P_AB"] #50:50 A:B
  means[,ncol(means)]         <- snp.param.indiv[,"MEAN_P_BB"] #homozygous B
  
  for (column in 2:(2+2*n.in.pools)) {
    
    if(column > 2 & column < (ncol(means)/2 + 1)) {
      means[,column] <- (snp.param.indiv[,"MEAN_P_AB"] * (column - 2) + 
                           snp.param.indiv[,"MEAN_P_AA"]  * ((ncol(means)/2 + 1) - column)) /
        (column - 2 + ncol(means)/2 + 1 - column)           
    }
    
    if(column > (ncol(means)/2 + 1) & column < ncol(means)) {
      means[,column] <- (snp.param.indiv[,"MEAN_P_BB"] * (column - (ncol(means)/2 + 1) ) + 
                           snp.param.indiv[,"MEAN_P_AB"]  * (ncol(means) - column)) /
        (column - (ncol(means)/2 + 1) + ncol(means) - column)
    }
  }
  
  #Generate standard deviations
  sds <- matrix(NA, nrow = nrow(snp.param.indiv), ncol = length(genotypes) + 1)
  sds <- as.data.frame(sds)
  colnames(sds) <- c("SNP_ID", paste("SD_P_",genotypes, sep=""))
  sds[,"SNP_ID"] <- snp.param.indiv[,"SNP_ID"]
  
  sds[,2]                 <- snp.param.indiv[,"SD_P_AA"] #homozygous A
  sds[,(ncol(sds)/2 + 1)] <- snp.param.indiv[,"SD_P_AB"] #50:50 A:B
  sds[,ncol(sds)]         <- snp.param.indiv[,"SD_P_BB"] #homozygous B
  
  for (column in 2:(2+2*n.in.pools)) {
    
    if(column > 2 & column < (ncol(sds)/2 + 1)) {
      sds[,column] <- sqrt((snp.param.indiv[,"SD_P_AB"]^2 * (column - 2) + 
                              snp.param.indiv[,"SD_P_AA"]^2  * ((ncol(sds)/2 + 1) - column)) /
                             (column - 2 + ncol(sds)/2 + 1 - column))      
    }
    
    if(column > (ncol(sds)/2 + 1) & column < ncol(sds)) {
      sds[,column] <- sqrt((snp.param.indiv[,"SD_P_BB"]^2 * (column - (ncol(sds)/2 + 1) ) + 
                              snp.param.indiv[,"SD_P_AB"]^2  * (ncol(sds) - column)) /
                             (column - (ncol(sds)/2 + 1) + ncol(sds) - column))
    }
  }
  
  snp.param.pools <- cbind(means, sds[,-1])
  snp.param.pools <- snp.param.pools[,order(c("A",genotypes, genotypes))] #reorder
  snp.param.pools[,"SNP_ID"] <- as.character(snp.param.pools[,"SNP_ID"])
  
  snp.param.pools$SNP_ID <- as.character(snp.param.pools$SNP_ID)
  # snp.param.pools <- merge(snp.param.pools, map[,c("SNP_ID", "A_ALLELE", "B_ALLELE")], by = "SNP_ID", all.x = TRUE)
  snp.param.pools$SNP_ID <- as.character(snp.param.pools$SNP_ID)
  alleles.by.snp$SNP_ID             <- as.character(alleles.by.snp$SNP_ID)
  snp.param.pools <- left_join(snp.param.pools, alleles.by.snp[,c("SNP_ID", "A_ALLELE", "B_ALLELE")], by = "SNP_ID")
  snp.param.pools <- snp.param.pools[order(snp.param.pools$SNP_ID),]
  
  return(snp.param.pools)
  
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export
 
genotypes.fun <- function(n) {
  genotypes <- matrix("X", nrow = 1+n, ncol = n)
  for (row in 1:(1+n)) {
    genotypes[row,] <- c(rep("A", n - row + 1), rep("B", (row - 1)))
  }
  genotypes <- apply(genotypes, 1, paste, collapse='')
  return(genotypes)
}
