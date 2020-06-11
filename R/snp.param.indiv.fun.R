# June 2020
# Matthew Hamilton

#' snp.param.indiv.fun
#' 
#' This function generates SNP parameters (specifically allelic proportion means, standard deviations and welch statistics) 
#' as described in the 'Estimation of SNP specific parameters' section of Henshall et al. (2014).  
#' If specified, it also generates signal intensity scatter plots and signal intensity against allelic proportion scatter plots 
#' (e.g. Figure 1 of Henshall et al. 2014) for each SNP.
#' 
#' @param snp.dat.indiv is a data frame with the following headings (class in parentheses):
#' \itemize{
#'  \item{'SAMPLE_ID' is the sample identifier.  Samples must be from diploid individuals (i.e. not pools) (integer).}
#'  \item{'SNP_ID' is the SNP identifier (character).} 
#'  \item{'INTENSITY_A' is the area/intensity for allele A (numeric).}
#'  \item{'INTENSITY_B' is the area/intensity for allele B (numeric).} 
#'  \item{'A_ALLELE' is the base represented by allele A (i.e. 'A', 'C', 'G' or 'T') (character).}
#'  \item{'B_ALLELE' is the base represented by allele B (i.e. 'A', 'C', 'G' or 'T') (character).}
#'  \item{'GENOTYPE' is the SNP genotype call (e.g. 'AT', 'TT').  NA if missing (character).}
#' } 

#' @param min.count is an integer (default = 2) representing the minimum number of each genotype (i.e. 'AA', 'AB', 'BB') for each SNP
#' required for the computation of SNP parameters.

#' @param min.intensity is a numeric variable (default = 0).  If the square root of the sum of INTENSITY_A squared and 
#' INTENSITY_B squared is less than min.intensity then this record is excluded in the computation of SNP parameters.
#' That is, observations that fall into an arc with a radius equal to min.intensity in the lower left of
#' signal intensity scatter plots are excluded.

#' @param gen.plots   is a logical variable.  If TRUE then signal intensity scatter plots and 
#' signal intensity against allelic proportion scatter plots for each SNP are generated and saved to disc 

#' @return 1. A data frame containing estimates of SNP sepecific parameters (refer to Henshall et al. 2014): 
#' \itemize{
#'  \item{'SNP_ID' is the SNP identifier.} 
#'  \item{'N_AA' is the count of homozygous A (AA) genotypes.}
#'  \item{'MEAN_P_AA' is the mean of allelic proportion for homozygous A genotypes.}
#'  \item{'SD_P_AA' is the standard deviation of allelic proportion for homozygous A genotypes.} 
#'  \item{'N_AB' is the count of heterozygous (AB) genotypes.}
#'  \item{'MEAN_P_AB' is the mean of allelic proportion for heterozygous (AB) genotypes.}
#'  \item{'SD_P_AB' is the standard deviation of allelic proportion for heterozygous (AB) genotypes.}  
#'  \item{'N_BB' is the count of homozygous B (BB) genotypes.}
#'  \item{'MEAN_P_BB' is the mean of allelic proportion for homozygous B genotypes.}
#'  \item{'SD_P_BB' is the standard deviation of allelic proportion for homozygous B genotypes.} 
#'  \item{'WELCH_A' is the welsh statistic for the interval between MEAN_P_AA and MEAN_P_AB.}
#'  \item{'WELCH_B' is the welsh statistic for the interval between MEAN_P_AB and MEAN_P_BB.}
#'  \item{'A_ALLELE_FREQ' is the A allele frequency computed from genotype counts.}
#'  \item{'B_ALLELE_FREQ' is the B allele frequency computed from genotype counts.}
#'  \item{'A_ALLELE' is the base represented by allele A (i.e. 'A', 'C', 'G' or 'T').}
#'  \item{'B_ALLELE' is the base represented by allele B (i.e. 'A', 'C', 'G' or 'T').}
#' } 

#' @return 2. '[SNP_ID].intensity.png' (if gen.plots = TRUE): signal intensity scatter plots saved in a sub-directory named 'Results/SNP_intensity_plots' 
#' in the working directory. 
#' @return 3. '[SNP_ID].prop.png' (if gen.plots = TRUE): signal intensity against allelic proportion scatter plots saved in a sub-directory named directory named
#' 'Results/SNP_allelic_prop_plots' in the working directory. See Figure 1 of Henshall et al. 2014

#' @examples
#' #Retrieve data for small worked example from Hamilton 2020
#' data(Ham.snp.dat.indiv)
#' 
#' #Compute SNP parameters
#' snp.param.indiv.fun(Ham.snp.dat.indiv)

#' @references Henshall JM, Dierens, L Sellars MJ (2014) Quantitative analysis of low-density SNP data for parentage assignment and estimation of family contributions to pooled samples. Genetics Selection Evolution 46, 51. https://doi 10.1186/s12711-014-0051-y 
#' @references Hamilton MG (2020) Maximum likelihood parentage assignment using quantitative genotypes

#' @export

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

snp.param.indiv.fun <- function(snp.dat.indiv, 
                          min.count = 2,             #snp.gen.param.fun
                          min.intensity = 0,
                          gen.plots = FALSE) { #pij.fun
  
  # load required packages
  if("dplyr" %in% installed.packages()[, "Package"] == FALSE) {install.packages("dplyr", repos='https://cran.csiro.au/')} 
  library(dplyr)
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
           "A_ALLELE", "B_ALLELE", "GENOTYPE") %in% colnames(snp.dat.indiv)) != 7) {
    stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B, A_ALLELE, B_ALLELE, GENOTYPE")
  }
  
  #Name columns and assign class
  #colnames(snp.dat.indiv) <- c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "A_ALLELE", "B_ALLELE", "GENOTYPE")
  snp.dat.indiv$SAMPLE_ID  <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$INTENSITY_A    <- as.numeric(snp.dat.indiv$INTENSITY_A)
  snp.dat.indiv$INTENSITY_B    <- as.numeric(snp.dat.indiv$INTENSITY_B)
  snp.dat.indiv$A_ALLELE   <- as.character(snp.dat.indiv$A_ALLELE)
  snp.dat.indiv$B_ALLELE   <- as.character(snp.dat.indiv$B_ALLELE)
  snp.dat.indiv$GENOTYPE  <- as.character(snp.dat.indiv$GENOTYPE)
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
                        "A_ALLELE", "B_ALLELE", "GENOTYPE")]
  
  if(sum(!is.integer(snp.dat.indiv$SAMPLE_ID)) > 0) {
    stop("SAMPLE_ID in snp.dat.indiv must be an integer.  Also check for missing values.")
  }
  if(sum(!is.character(snp.dat.indiv$SNP_ID)) > 0) {
    stop("SNP_ID in snp.dat.indiv must be a character.  Check for missing values.")
  }
  if(sum(!is.numeric(snp.dat.indiv$INTENSITY_A)) > 0) {
    stop("INTENSITY_A in snp.dat.indiv must be numeric.  Also check for missing values.")
  }
  if(sum(!is.numeric(snp.dat.indiv$INTENSITY_B)) > 0) {
    stop("INTENSITY_B in snp.dat.indiv must be numeric.  Also check for missing values.")
  }
  if(sum(!is.character(snp.dat.indiv$A_ALLELE)) > 0) {
    stop("A_ALLELE in snp.dat.indiv must be a character.  Check for missing values.")
  }
  if(sum(!is.character(snp.dat.indiv$B_ALLELE)) > 0) {
    stop("B_ALLELE in snp.dat.indiv must be a character.  Check for missing values.")
  }
  if(sum(!is.character(snp.dat.indiv$GENOTYPE)) > 0) {
    stop("GENOTYPE in snp.dat.indiv must be a character.  Check for missing values.")
  }
  
  #Check for duplicated records in snp.dat.indiv
  indiv.snp <- paste(snp.dat.indiv$SAMPLE_ID,snp.dat.indiv$SNP_ID, sep=".")
  if(sum(duplicated(indiv.snp)) > 0) {
    stop("SAMPLE_ID and SNP_ID combinations are not unique in snp.dat.indiv.  Delete duplicates or recode SAMPLE_ID.")
  }
  rm(indiv.snp)
  
  #Check that there is only one nucleotide in column A_ALLELE for each SNP
  tmp.1 <-  unique(snp.dat.indiv[,c("SNP_ID","A_ALLELE")])
  if(sum(unique(snp.dat.indiv[,"SNP_ID"]) != tmp.1[,c("SNP_ID")]) > 0) {
    stop("There is more than one nucleotide in column A_ALLELE for at least one SNP")
  }
  rm(tmp.1)
  
  #Check that there is only one nucleotide in column B_ALLELE for each SNP  
  tmp.1 <-  unique(snp.dat.indiv[,c("SNP_ID","B_ALLELE")])
  if(sum(unique(snp.dat.indiv[,"SNP_ID"]) != tmp.1[,c("SNP_ID")]) > 0) {
    stop("There is more than one nucleotide in column B_ALLELE for at least one SNP")
  }
  rm(tmp.1)
  
  #Check that A_ALLELE and B_ALLELE is one of A, C, G or T
  if(!sum(unique(snp.dat.indiv[,"A_ALLELE"]) %in% c("A", "C", "G", "T")) == length(unique(snp.dat.indiv[,"A_ALLELE"]))){
    stop("A_ALLELE must be one of A, C, G or T.  Check for missing values")
  }
  if(!sum(unique(snp.dat.indiv[,"B_ALLELE"]) %in% c("A", "C", "G", "T")) == length(unique(snp.dat.indiv[,"B_ALLELE"]))){
    stop("B_ALLELE must be one of A, C, G or T.  Check for missing values")
  }
  
  #Check that GENOTYPE is comprised of 2 characters
  if(!sum(is.na(snp.dat.indiv[,"GENOTYPE"]) | nchar(snp.dat.indiv[,"GENOTYPE"]) == 2) == nrow(snp.dat.indiv)) {
    stop("If not missing, GENOTYPE must be comprised of 2 characters")
  }
  
  #Check that GENOTYPE is comprised of A_ALLELE and/or B_ALLELE
  if(!sum(is.na(snp.dat.indiv[,"GENOTYPE"]) |
          ((substring(snp.dat.indiv[,"GENOTYPE"],1,1) == snp.dat.indiv[,"A_ALLELE"] |
            substring(snp.dat.indiv[,"GENOTYPE"],1,1) == snp.dat.indiv[,"B_ALLELE"] ) &
           (substring(snp.dat.indiv[,"GENOTYPE"],2,2) == snp.dat.indiv[,"A_ALLELE"] |
            substring(snp.dat.indiv[,"GENOTYPE"],2,2) == snp.dat.indiv[,"B_ALLELE"] ))) ==
     nrow(snp.dat.indiv)) {
    stop("If not missing, GENOTYPE must be comprised of nucleotides in A_ALLELE and B_ALLELE for each SNP")
  }
  
  #End Checks####################################################
  wd <- getwd()
  
  if(gen.plots) {
    
    dir.create(file.path(wd, "SNP_outputs"), showWarnings = FALSE)
    
    setwd(file.path(wd, "SNP_outputs"))
    
    snp.intensity.plot.fun(snp.dat.indiv = snp.dat.indiv) 
  }
  
  snp.param.indiv <- snp.gen.param.fun(snp.dat.indiv = snp.dat.indiv,
                                 min.count = min.count,
                                 min.intensity = min.intensity)
  setwd(wd)
  
  return(snp.param.indiv)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

#Define snp.intensity.plot.fun function

snp.intensity.plot.fun <- function(snp.dat.indiv, 
                                   min.intensity = 0) {
  
  #Generates SNP intensity scatter plots from genotype calls intensity vs allelic proportion
  #plots akin to Figure 1 of Henshall et al. 2014
  
  #Args##########################################
  # snp.dat.indiv: Data frame
  #             1. SAMPLE_ID is the individual identifier
  #             2. SNP_ID   is the SNP identifier
  #             3. INTENSITY_A   is the area/intensity for allele A
  #             4. INTENSITY_B   is the area/intensity for allele B
  #             7. A_ALLELE  is the base represented by allele A
  #             8. B_ALLELE  is the base represented by allele B
  #             9. GENOTYPE is the SNP genotype call
  
  # min.intensity      Number used in pij.fun. If sqrt((snp.dat.indiv$INTENSITY_A)^2 +
  #              (snp.dat.indiv$INTENSITY_B)^2) less than this value
  #              then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  #              Essentially removes observations that fall into the lower left of INTENSITY_A
  #              by INTENSITY_B scatter plot.
  
  #Returns##########################################
  # [SNP_ID].intensity.png:  Plot of intensity A vs intensity B returned in a directory named "SNP_intensity_plots"
  #                            in the working directory. 
  # [SNP_ID].prop.png:  Plot of allelic proportion vs intensity returned in a directory named 
  #                     "SNP_allelic_prop_plots" in the working directory. See Figure 1 of Henshall et al. 2014
  
  print("Running snp.intensity.plot.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
           "A_ALLELE", "B_ALLELE", "GENOTYPE") %in% colnames(snp.dat.indiv)) != 7) {
    stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B, A_ALLELE, B_ALLELE, GENOTYPE")
  }
  
  #Name columns and assign class
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
                        "A_ALLELE", "B_ALLELE", "GENOTYPE")]
  snp.dat.indiv$SAMPLE_ID <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$INTENSITY_A <- as.numeric(snp.dat.indiv$INTENSITY_A)
  snp.dat.indiv$INTENSITY_B <- as.numeric(snp.dat.indiv$INTENSITY_B)
  snp.dat.indiv$A_ALLELE <- as.character(snp.dat.indiv$A_ALLELE)
  snp.dat.indiv$B_ALLELE <- as.character(snp.dat.indiv$B_ALLELE)
  snp.dat.indiv$GENOTYPE <- as.character(snp.dat.indiv$GENOTYPE)  
  
  #Check for duplicated records in snp.dat.indiv
  indiv.snp <- paste(snp.dat.indiv$SAMPLE_ID,snp.dat.indiv$SNP_ID, sep=".")
  if(sum(duplicated(indiv.snp)) > 0) {
    stop("SAMPLE_ID and SNP_ID combinations are not unique in snp.dat.indiv.  Delete duplicates or recode SAMPLE_ID.")
  }
  rm(indiv.snp)
  
  #Get allelic proportion
  snp.allelic.prop <- pij.fun(snp.dat.indiv = snp.dat.indiv, min.intensity = min.intensity)
  print("Still running snp.intensity.plot.fun")
  
  #Rename columns
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "SAMPLE_ID"]        <- "SAMPLE_ID"
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "ALLELIC_PROP"]     <- "ALLELIC_PROP_INDIV"  
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "INTENSITY"]        <- "INTENSITY_INDIV"  
  
  # snp.dat.indiv <- merge(snp.dat.indiv, snp.allelic.prop, by = c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B"), all.x = TRUE)
  snp.dat.indiv <- left_join(snp.dat.indiv, snp.allelic.prop, by = c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B"))
  
  #Generate directory for plots if does not exist
  if (!file.exists("SNP_intensity_plots")){
    dir.create(file.path(getwd(), "SNP_intensity_plots"))
  }
  
  if (!file.exists("SNP_allelic_prop_plots")){
    dir.create(file.path(getwd(), "SNP_allelic_prop_plots"))
  }
  
  setwd(file.path(getwd(), "SNP_intensity_plots"))
  
  #Plot parameters
  plot.width <-  680
  plot.height <-  680
  plot.pointsize <- 24
  
  #Intensity plots  
  
  #Generate plots for each SNP
  for (snp in unique(snp.dat.indiv$SNP_ID)) { 
    
    indiv.snp.dat.indiv <- snp.dat.indiv[snp.dat.indiv$SNP_ID == snp,]
    
    png(filename = paste(snp,".intensity.png",sep=""), 
        width = 1*plot.width, 
        height = 1*plot.height,
        pointsize = plot.pointsize)
    
    par(mfrow=c(1,1)) #1 plot
    
    plot.indiv.fun(plot.type      = "Intensity",
                   indiv.snp.dat.indiv  = indiv.snp.dat.indiv)
    
    dev.off()
  }
  
  #Allelic proportion plots
  
  setwd("../SNP_allelic_prop_plots")
  
  #Generate plots for each SNP
  for (snp in unique(snp.dat.indiv$SNP_ID)) { 
    indiv.snp.dat.indiv <- snp.dat.indiv[snp.dat.indiv$SNP_ID == snp,]
    
    png(filename = paste(snp,".prop.png",sep=""), 
        width = 1*plot.width, 
        height = 1*plot.height,
        pointsize = plot.pointsize)
    
    par(mfrow=c(1,1)) #1 plot
    
    plot.indiv.fun(plot.type      = "Prop", 
                   indiv.snp.dat.indiv  = indiv.snp.dat.indiv)
    
    dev.off()
  }
  
  setwd("../")
  
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

snp.gen.param.fun <- function(snp.dat.indiv, 
                              min.count = 2,
                              min.intensity = 0) {
  #Description:
  # Generates SNP parameters
  
  #Args##########################################
  # snp.dat.indiv: Data frame
  #              SAMPLE_ID is the individual identifier
  #              SNP_ID   is the SNP identifier
  #              INTENSITY_A   is the area/intensity for allele A
  #              INTENSITY_B   is the area/intensity for allele B
  #              A_ALLELE  is the base represented by allele A
  #              B_ALLELE  is the base represented by allele B
  #              GENOTYPE is the SNP genotype call
  
  # min.count:   Number.  Minimum number of each genotype (AA, AB, BB) for each SNP required
  #              for the calculation of allelic proportion means, standard deviations and 
  #              welch statistics.
  
  # min.intensity      Number used in pij.fun. If sqrt((snp.dat.indiv$INTENSITY_A)^2 +
  #              (snp.dat.indiv$INTENSITY_B)^2) less than this value
  #              then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  #              Removes observations that fall into the lower left of INTENSITY_A
  #              by INTENSITY_B scatter plot.
  
  #Returns##########################################
  
  #   Data frame. See "Estimation of SNP sepecific parameters" page 3 of Henshall et al 2014 
  #              SNP_ID        is the SNP identifier, 
  #              N_AA          is the count homozygotes for allele A
  #              MEAN_P_AA     is the mean of allelic proportion (homozygous allele A), 
  #              SD_P_AA       is the standard deviation of allelic proportion (homozygous allele A), 
  #              N_AB          is the count heterozygotes
  #              MEAN_P_AB     is the mean of allelic proportion (heterozygous), 
  #              SD_P_AB       is the standard deviation of allelic proportion (heterozygous), 
  #              N_BB          is the count homozygotes for allele B
  #              MEAN_P_BB     is the mean of allelic proportion (homozygous allele B), 
  #              SD_P_BB       is the standard deviation of allelic proportion (homozygous allele B),
  #              WELCH_A       is the welsh statistic for the interval between AA and AB
  #              WELCH_B       is the welsh statistic for the interval between AB and BB
  #              A_ALLELE_FREQ is the A allele frequency derived from counts
  #              B_ALLELE_FREQ is the B allele frequency derived from counts
  #              A_ALLELE  is the base represented by allele A
  #              B_ALLELE  is the base represented by allele B
  
  print("Running snp.gen.param.fun")
  
  #Get allelic proportion
  snp.allelic.prop <- pij.fun(snp.dat.indiv = snp.dat.indiv, min.intensity = min.intensity)
  print("Still running snp.gen.param.fun")
  
  #Rename columns
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "SAMPLE_ID"]               <- "SAMPLE_ID"
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "ALLELIC_PROP"]     <- "ALLELIC_PROP_INDIV"  
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "INTENSITY"]        <- "INTENSITY_INDIV"  
  
  if(identical(snp.dat.indiv$SAMPLE_ID,snp.allelic.prop$SAMPLE_ID) &
     identical(snp.dat.indiv$SNP_ID,snp.allelic.prop$SNP_ID)) {
    snp.dat.indiv <- cbind(snp.dat.indiv, snp.allelic.prop[,-c(1:2)])
  } else {
    stop("SAMPLE_ID and SNP_ID columns of pij.fun output do not match those of snp.dat.indiv.  Not sure why.")
  }
  
  #Change missing values
  snp.dat.indiv[is.na(snp.dat.indiv[,"GENOTYPE"]), "GENOTYPE"] <- "No.call"
  
  snp.dat.indiv$ALLELES <- NA
  #Generate column of allele genotypes AA, AB, BB
  
  snp.dat.indiv[snp.dat.indiv[,"GENOTYPE"] == 
            paste(snp.dat.indiv[,"A_ALLELE"], snp.dat.indiv[,"A_ALLELE"], sep = ""),"ALLELES"] <- "AA"  
  snp.dat.indiv[snp.dat.indiv[,"GENOTYPE"] == 
            paste(snp.dat.indiv[,"B_ALLELE"], snp.dat.indiv[,"B_ALLELE"], sep = ""),"ALLELES"] <- "BB"  
  snp.dat.indiv[snp.dat.indiv[,"GENOTYPE"] == 
            paste(snp.dat.indiv[,"A_ALLELE"], snp.dat.indiv[,"B_ALLELE"], sep = ""),"ALLELES"] <- "AB"
  snp.dat.indiv[snp.dat.indiv[,"GENOTYPE"] == 
            paste(snp.dat.indiv[,"B_ALLELE"], snp.dat.indiv[,"A_ALLELE"], sep = ""),"ALLELES"] <- "AB"
  
  #Change missing values back to NA
  snp.dat.indiv[snp.dat.indiv[,"GENOTYPE"] == "No.call", "GENOTYPE"] <- NA
  
  #reorder
  snp.dat.indiv <- snp.dat.indiv[order(snp.dat.indiv$SAMPLE_ID, decreasing = FALSE), ] 
  
  #Get A and B alleles
  alleles <- unique(snp.dat.indiv[,c("SNP_ID", "A_ALLELE", "B_ALLELE")])
  
  #Loop through SNP
  snp.param.indiv <- NULL
  for(snp in unique(snp.dat.indiv$SNP_ID)) {
    
    #Get current SNP data
    current.snp.dat.indiv <- snp.dat.indiv[snp.dat.indiv[,"SNP_ID"] == snp,]
    
    #Get data for each combination of alleles
    
    #AA
    p.aa <- current.snp.dat.indiv[current.snp.dat.indiv[,"ALLELES"] == "AA","ALLELIC_PROP_INDIV"]
    p.aa <- p.aa[!is.na(p.aa)] #remove if NA
    if(length(p.aa) == 0) {p.aa <- NULL}
    
    #AB
    p.ab <- current.snp.dat.indiv[current.snp.dat.indiv[,"ALLELES"] == "AB","ALLELIC_PROP_INDIV"]
    p.ab <- p.ab[!is.na(p.ab)] #remove if NA
    if(length(p.ab) == 0) {p.ab <- NULL}
    
    #BB
    p.bb <- current.snp.dat.indiv[current.snp.dat.indiv[,"ALLELES"] == "BB","ALLELIC_PROP_INDIV"]
    p.bb <- p.bb[!is.na(p.bb)] #remove if NA
    if(length(p.bb) == 0) {p.bb <- NULL}
    
    #Compute allelic proportion counts, means and standard deviations for current SNP
    
    #Counts
    count.aa <- length(p.aa)
    count.ab <- length(p.ab)
    count.bb <- length(p.bb)
    
    #AA
    if(count.aa >= min.count & count.aa >= 1) {  
      mean.p.aa <- mean(p.aa)
    } else {
      mean.p.aa <- NA
    }
    
    if(count.aa >= min.count & count.aa >= 2) {  
      sd.p.aa   <- sd(p.aa)
    } else {
      sd.p.aa <- NA
    }
    
    #AB
    if(count.ab >= min.count & count.ab >= 1) {  
      mean.p.ab <- mean(p.ab)
    } else {
      mean.p.ab <- NA
    }
    
    if(count.ab >= min.count & count.ab >= 2) {  
      sd.p.ab   <- sd(p.ab)
    } else {
      sd.p.ab <- NA
    }
    
    #BB
    if(count.bb >= min.count & count.bb >= 1) {  
      mean.p.bb <- mean(p.bb)
    } else {
      mean.p.bb <- NA
    }
    
    if(count.bb >= min.count & count.bb >= 2) {  
      sd.p.bb   <- sd(p.bb)
    } else {
      sd.p.bb <- NA
    }
    
    #Calculate welch statistics for current SNP
    
    #WELCH_A
    if(count.aa >= min.count & count.ab >= min.count &
       count.aa >= 2 & count.ab >= 2) {
      welch.a   <- t.test(p.ab, p.aa)$statistic
    } else {
      welch.a <- NA      
    }
    
    #WELCH_B
    if(count.bb >= min.count & count.ab >= min.count &
       count.bb >= 2 & count.ab >= 2) {
      welch.b   <- t.test(p.bb, p.ab)$statistic   
    } else {
      welch.b <- NA
    }
    
    #A_ALLELE_FREQ
    a.allele.freq <- (count.aa + count.ab/2) / (count.aa + count.ab + count.bb)
    #B_ALLELE_FREQ
    b.allele.freq <- (count.bb + count.ab/2) / (count.aa + count.ab + count.bb)
    
    #Generate data frame (one row) with current SNP parameters
    current.snp.dat.indiv <- data.frame(
      SNP_ID    = snp,
      N_AA      = count.aa,
      MEAN_P_AA = mean.p.aa,
      SD_P_AA   = sd.p.aa,
      N_AB      = count.ab,
      MEAN_P_AB = mean.p.ab,
      SD_P_AB   = sd.p.ab,
      N_BB      = count.bb,
      MEAN_P_BB = mean.p.bb,
      SD_P_BB   = sd.p.bb,  
      WELCH_A   = welch.a,
      WELCH_B   = welch.b,
      A_ALLELE_FREQ = a.allele.freq,
      B_ALLELE_FREQ = b.allele.freq
    )
    snp.param.indiv <- rbind(snp.param.indiv,current.snp.dat.indiv)
    
    rm(p.aa, p.ab, p.bb, current.snp.dat.indiv, mean.p.aa, sd.p.aa, mean.p.ab,
       sd.p.ab, mean.p.bb, sd.p.bb, welch.a, welch.b, 
       a.allele.freq, b.allele.freq)
  }
  
  colnames(snp.param.indiv) <- c("SNP_ID", "N_AA", "MEAN_P_AA", "SD_P_AA", 
                           "N_AB", "MEAN_P_AB", "SD_P_AB", 
                           "N_BB", "MEAN_P_BB", "SD_P_BB", 
                           "WELCH_A", "WELCH_B", 
                           "A_ALLELE_FREQ", "B_ALLELE_FREQ")
  
  #Add A_ALLELE and B_ALLELE columns
  #snp.param.indiv <- merge(snp.param.indiv, alleles, by = "SNP_ID", all.x = TRUE)
  alleles$SNP_ID    <- as.character(alleles$SNP_ID)
  snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  snp.param.indiv <- left_join(snp.param.indiv, alleles, by = "SNP_ID")
  
  snp.param.indiv$SNP_ID <- as.character(snp.param.indiv$SNP_ID)
  snp.param.indiv <- snp.param.indiv[order(snp.param.indiv$SNP_ID),]
  
  return(snp.param.indiv)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

pij.fun <- function(snp.dat.indiv, min.intensity = 0, min.intensity.old = FALSE) {
  
  #Equation 1 of Henshall et al. 2014
  
  #Args##########################################
  # snp.dat.indiv: Data frame  
  #            1. SAMPLE_ID is the individual (or pool) identifier, 
  #            2. SNP_ID is the SNP identifier, 
  #            3. INTENSITY_A is the area/intensity for allele A
  #            4. INTENSITY_B is the area/intensity for allele B 
  
  # min.intensity      Number. If sqrt((snp.dat.indiv$INTENSITY_A)^2 +
  #              (snp.dat.indiv$INTENSITY_B)^2) less than this value
  #              then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  #              Essentially removes observations that fall into the lower left of INTENSITY_A
  #              by INTENSITY_B scatter plot.
  
  # min.intensity.old: Logical.  Generally FALSE.  If true the the rule applied to min.intensity is changed to
  #                    If (snp.dat.indiv$INTENSITY_A + snp.dat.indiv$INTENSITY_B) less than min.intensity
  #                    then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  
  #Returns##########################################
  
  #snp.allelic.prop: Data frame.  Returns allelic proportion for each SNP
  #            1. SAMPLE_ID is the identifier (pool or individual)
  #            2. SNP_ID is the SNP identifier
  #            3. INTENSITY_A is the area/intensity for allele A
  #            4. INTENSITY_B is the area/intensity for allele B 
  #            5. ALLELIC_PROP is the allelic proportion
  #            6. INTENSITY is the hypotenuse: sqrt(INTENSITY_A^2 + INTENSITY_B^2)
  
  #Name columns and assign class
  #colnames(snp.dat.indiv) <- c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")
  
  print("Running pij.fun")
  
  colnames(snp.dat.indiv)[colnames(snp.dat.indiv) == "SAMPLE_ID"] <- "SAMPLE_ID"
  colnames(snp.dat.indiv)[colnames(snp.dat.indiv) == "SAMPLE_ID"]  <- "SAMPLE_ID"
  
  snp.dat.indiv$SAMPLE_ID        <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$INTENSITY_A    <- as.numeric(snp.dat.indiv$INTENSITY_A)
  snp.dat.indiv$INTENSITY_B    <- as.numeric(snp.dat.indiv$INTENSITY_B)
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")]
  
  # Set allelic proportion to missing according to value of min.intensity
  if(min.intensity.old == FALSE) {
    snp.dat.indiv[sqrt(snp.dat.indiv[,"INTENSITY_A"]^2 + snp.dat.indiv[,"INTENSITY_B"]^2) < min.intensity &
              !is.na(snp.dat.indiv[,"INTENSITY_A"]^2 + snp.dat.indiv[,"INTENSITY_B"]^2) , 
            c("INTENSITY_A","INTENSITY_B")] <- c(NA, NA)
  } else {
    snp.dat.indiv[(snp.dat.indiv[,"INTENSITY_A"] + snp.dat.indiv[,"INTENSITY_B"]) < min.intensity, 
            c("INTENSITY_A","INTENSITY_B")] <- c(NA, NA)    
  }
  
  #Equation 1 of Henshall et al. 2014
  snp.dat.indiv$ALLELIC_PROP <- atan(snp.dat.indiv$INTENSITY_B / snp.dat.indiv$INTENSITY_A)/(pi/2) 
  snp.dat.indiv[snp.dat.indiv[,"INTENSITY_A"] < 0 & !is.na(snp.dat.indiv[,"INTENSITY_A"]),"ALLELIC_PROP"] <- 
    2+snp.dat.indiv[snp.dat.indiv[,"INTENSITY_A"] < 0 & !is.na(snp.dat.indiv[,"INTENSITY_A"]),"ALLELIC_PROP"] #NOTE where INTENSITY_A < 0 ALLELIC_PROP is negative - need to add 2
  snp.dat.indiv$INTENSITY <- sqrt(snp.dat.indiv$INTENSITY_B^2 + snp.dat.indiv$INTENSITY_A^2)
  
  snp.allelic.prop <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
                                 "ALLELIC_PROP", "INTENSITY")]
  
  return(snp.allelic.prop)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

plot.indiv.fun <- function(plot.type, indiv.snp.dat.indiv) {
  
  #Generates single intensity plot.  Required by snp.intensity.plot.fun.
  
  #Args##########################################
  # plot.type:  Character string.  "Intensity" or "Prop" If "Intensity" then 
  #             a plot of Intensity.A by Intensity B is generated.  
  #             If "Prop" then plot of intensity against allelic proportion. 
  # indiv.snp.dat.indiv: Data frame
  #              1. SAMPLE_ID is the individual identifier
  #              2. SNP_ID   is the SNP identifier
  #              3. INTENSITY_A   is the area/intensity for allele A
  #              4. INTENSITY_B   is the area/intensity for allele B
  #              7. A_ALLELE  is the base represented by allele A
  #              8. B_ALLELE  is the base represented by allele B
  #              9. GENOTYPE is the SNP genotype call
  #              12. ALLELIC_PROP_INDIV is the allelic proportion (adjusted for uncertainty; See Equation 1 of Henshall et al. 2014)   
  #              13. INTENSITY_INDIV is the intensity not adjusted for uncertainty (see Figure 1 of Henshall et al 2014)      
  
  #Returns##########################################
  # scatterplot:  Plot of intensity A vs intensity B returned 
  #               in a directory named "snp.scatter" in the working directory. 
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
           "A_ALLELE", "B_ALLELE", "GENOTYPE", "ALLELIC_PROP_INDIV", 
           "INTENSITY_INDIV") %in% colnames(indiv.snp.dat.indiv)) != 9) {
    stop("indiv.snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B, A_ALLELE, B_ALLELE, GENOTYPE, ALLELIC_PROP_INDIV, INTENSITY_INDIV")
  }
  
  #Name columns and assign class
  #colnames(indiv.snp.dat.indiv)  <- c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "A_ALLELE", "B_ALLELE", "GENOTYPE")
  indiv.snp.dat.indiv$SAMPLE_ID   <- as.integer(indiv.snp.dat.indiv$SAMPLE_ID)
  indiv.snp.dat.indiv$SNP_ID      <- as.character(indiv.snp.dat.indiv$SNP_ID)
  indiv.snp.dat.indiv$INTENSITY_A <- as.numeric(indiv.snp.dat.indiv$INTENSITY_A)
  indiv.snp.dat.indiv$INTENSITY_B <- as.numeric(indiv.snp.dat.indiv$INTENSITY_B)
  indiv.snp.dat.indiv$A_ALLELE    <- as.character(indiv.snp.dat.indiv$A_ALLELE)
  indiv.snp.dat.indiv$B_ALLELE    <- as.character(indiv.snp.dat.indiv$B_ALLELE)
  indiv.snp.dat.indiv$GENOTYPE    <- as.character(indiv.snp.dat.indiv$GENOTYPE)           
  indiv.snp.dat.indiv$ALLELIC_PROP_INDIV  <- as.numeric(indiv.snp.dat.indiv$ALLELIC_PROP_INDIV) 
  indiv.snp.dat.indiv$INTENSITY_INDIV     <- as.numeric(indiv.snp.dat.indiv$INTENSITY_INDIV)    
  indiv.snp.dat.indiv <- indiv.snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
                                    "A_ALLELE", "B_ALLELE", "GENOTYPE", 
                                    "ALLELIC_PROP_INDIV", 
                                    "INTENSITY_INDIV")]
  
  #Check that plot.type = "Intensity",  "Prop"
  if(!(plot.type %in% c("Intensity", "Prop"))) {
    stop("plot.type must equal Intensity or Prop when calling the snp.intensity.plot.fun function")
  }
  
  #identify SNP
  snp <- unique(indiv.snp.dat.indiv$SNP_ID)
  
  #Define plot variables and data
  if(plot.type == "Intensity") {
    indiv.snp.dat.indiv$X.AXIS <- indiv.snp.dat.indiv$INTENSITY_A
    x.title <- "Intensity A"
    indiv.snp.dat.indiv$Y.AXIS <- indiv.snp.dat.indiv$INTENSITY_B    
    y.title <- "Intensity B"
    main.title <- paste(snp)
    lim.x <- c(0,max(c(indiv.snp.dat.indiv$X.AXIS, indiv.snp.dat.indiv$Y.AXIS), na.rm = TRUE))
    lim.y <- lim.x
  }
  
  if(plot.type == "Prop") {
    indiv.snp.dat.indiv$X.AXIS <- indiv.snp.dat.indiv$ALLELIC_PROP_INDIV
    x.title <- "Allelic proportion"
    indiv.snp.dat.indiv$Y.AXIS <- indiv.snp.dat.indiv$INTENSITY_INDIV    
    y.title <- "Intensity"
    main.title <- paste(snp)
    lim.x <- c(0,1)
    lim.y <- c(0,1.35*max(indiv.snp.dat.indiv$Y.AXIS, na.rm = TRUE)) #allow room for legend
  }
  
  #define x and y axis limits
  
  #Make Genotype a factor and recategorise NAs
  
  indiv.snp.dat.indiv[,"GENOTYPE"] <- as.character(indiv.snp.dat.indiv[,"GENOTYPE"])
  #Change NA to '-' for sorting
  indiv.snp.dat.indiv[is.na(indiv.snp.dat.indiv[,"GENOTYPE"]),"GENOTYPE"] <- "-"
  #Reorder: "-" always last
  indiv.snp.dat.indiv <- indiv.snp.dat.indiv[order(indiv.snp.dat.indiv$GENOTYPE, decreasing = TRUE), ] 
  #Reorder: heterozygote always first
  indiv.snp.dat.indiv <- indiv.snp.dat.indiv[order(nchar(as.vector(indiv.snp.dat.indiv[,"GENOTYPE"])), 
                                       decreasing = TRUE), ]  
  
  #Change genotype names
  indiv.snp.dat.indiv[indiv.snp.dat.indiv[,"GENOTYPE"] == "-","GENOTYPE"] <- "No call"
  indiv.snp.dat.indiv[indiv.snp.dat.indiv[,"GENOTYPE"] == "A","GENOTYPE"] <- "AA"
  indiv.snp.dat.indiv[indiv.snp.dat.indiv[,"GENOTYPE"] == "C","GENOTYPE"] <- "CC"
  indiv.snp.dat.indiv[indiv.snp.dat.indiv[,"GENOTYPE"] == "G","GENOTYPE"] <- "GG"
  indiv.snp.dat.indiv[indiv.snp.dat.indiv[,"GENOTYPE"] == "T","GENOTYPE"] <- "TT"
  
  #As factor (levels in order)
  indiv.snp.dat.indiv$GENOTYPE <- factor(indiv.snp.dat.indiv$GENOTYPE, 
                                   levels = unique(indiv.snp.dat.indiv$GENOTYPE))
  
  #Define colours depending on the levels of GENOTYPE
  #Colourblind friendly colour pallet 
  #From http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/
  #additional colours: "#E69F00" #Orange;         "#56B4E9" #Sky blue; 
  #                    "#009E73" #Bluish green;   "#F0E442" #Yellow
  #                    "#CC79A7",#Redish purple ; "#D55E00", #vermillion ; 
  #                    "#0072B2", #Blue  ;        "#000000" #Black
  
  colours <- NA
  #two homozygotes, one heterozygote, one missing
  if(length(levels(indiv.snp.dat.indiv$GENOTYPE)) == 4) {
    colours <- c("#009E73", #Bluish green;
                 "#E69F00", #Orange;
                 "#CC79A7",#Redish purple 
                 "#000000") #Black
  }
  
  #two homozygotes, one heterozygote
  if(length(levels(indiv.snp.dat.indiv$GENOTYPE)) == 3) {
    if(!("No call" %in% levels(indiv.snp.dat.indiv$GENOTYPE))) { #none missing
      colours <- c("#009E73", #Bluish green;
                   "#E69F00", #Orange;
                   "#CC79A7")#Redish purple
    }
  }
  
  #two (assume one hetero and one homo) and one missing 
  if(length(levels(indiv.snp.dat.indiv$GENOTYPE)) == 3) {
    if(("No call" %in% levels(indiv.snp.dat.indiv$GENOTYPE))) { #missing
      colours <- c("#009E73", #Bluish green;
                   "#E69F00", #Orange; 
                   "#000000") #Black 
    }
  }
  
  #two (assume one hetero and one homo) 
  if(length(levels(indiv.snp.dat.indiv$GENOTYPE)) == 2) {
    if(!("No call" %in% levels(indiv.snp.dat.indiv$GENOTYPE))) { #none missing
      colours <- c("#009E73", #Bluish green;
                   "#E69F00") #Orange;
    }
  }
  
  #one (asssume homo) and one missing 
  if(length(levels(indiv.snp.dat.indiv$GENOTYPE)) == 2) {
    if(("No call" %in% levels(indiv.snp.dat.indiv$GENOTYPE))) { #missing
      colours <- c("#009E73", #Bluish green;
                   "#000000") #Black 
    }
  }
  
  #one homozygote
  if(length(levels(indiv.snp.dat.indiv$GENOTYPE)) == 1) {
    colours <- c("#009E73") #Bluish green;
  }
  
  #Generate plot
  plot(indiv.snp.dat.indiv$X.AXIS,indiv.snp.dat.indiv$Y.AXIS,
       type="p",pch = 16, cex = 0.5,
       xlim = lim.x,
       ylim = lim.y,
       main= main.title,
       xlab = x.title,
       ylab = y.title,
       col = colours[indiv.snp.dat.indiv$GENOTYPE])
  
  #Get legend text prior to plotting
  if(plot.type == "Intensity") {
    legend.txt <- levels(indiv.snp.dat.indiv$GENOTYPE)
  }
  
  
  if(plot.type == "Prop") {
    
    if(plot.type == "Prop") {
      tmp.dat.all <- indiv.snp.dat.indiv[, "ALLELIC_PROP_INDIV"]
    } 
    
    legend.txt <- NULL
    for(geno in levels(indiv.snp.dat.indiv$GENOTYPE)) {
      tmp.dat    <- tmp.dat.all[indiv.snp.dat.indiv[,"GENOTYPE"] == geno]
      x.mean     <- round(mean(tmp.dat, na.rm = TRUE),2)
      x.sd       <- round(sd(tmp.dat, na.rm = TRUE), 2)
      tmp.legend <- (paste(geno," (mu = ",x.mean, ", sd = ",x.sd,")",sep=""))
      legend.txt <- c(legend.txt,tmp.legend)
      
    }
  }
  
  legend(x="topright",
         legend = legend.txt, 
         fill=colours,
         cex=0.8,
         bty = "n")
}


