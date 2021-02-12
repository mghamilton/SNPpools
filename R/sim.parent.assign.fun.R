# Feb 2021
# Matthew Hamilton

#' sim.parent.assign.fun
#'
#' @description 
#' This function adopts a stochastic simulation approach to determine the proportion of correct assignments and, for 
#' maximum likelihood approaches, the critical delta LOD values.  For each repetition, 
#' 'snp.dat.indiv', 'snp.dat.pools', 'snp.param.indiv', 'snp.param.pools', 'fam.set.combns' and 'fam.set.combns.by.pool' data frames 
#' for one pooled DNA sample are generated from user-defined 'ped', 'map', 'true.snp.param.indiv' and 'sim.fam.sets' data frames. 
#' Parentage is assigned for each simulated pool using the parent.assign.fun.

#' @param n_repetitions is a integer variable defining the number of repetitions in the simulation

#' @param ped is a data frame and a conventional pedigree file with one additional column.  
#' It must include all SAMPLE_IDs used to construct snp.param.indiv and all possible parents of pooled samples. 
#' It contains the following headings (class in parentheses):
#' \itemize{
#'  \item{'SAMPLE_ID' is the individual (i.e. not a pooled sample) sample identifier.  Individuals 
#'  with no true SAMPLE_ID should be assigned a dummy SAMPLE_ID.(integer).} 
#'  \item{'SIRE_ID' is the SAMPLE_ID of the sire (0 if unknown) (integer).} 
#'  \item{'DAM_ID' is the SAMPLE_ID of the dam  (0 if unknown) (integer).} 
#'  \item{'SAMPLED' if TRUE individual used to generate snp.param.indiv (logical).} 
#' }
#' @param map is a data frame and genetic map identifying the position of SNP.  
#' \itemize{
#'  \item{'CHROMOSOME' is the chromosome number.  To assume that SNP are not linked provided a unique CHROMOSOME number
#'   for each SNP_ID (integer).} 
#'  \item{'SNP_ID' is the SNP identifier (ordered by physical position within chromosome) (character).} 
#'  \item{'GENETIC_POSITION' is the SNP genetic position in Morgans (numeric). To assume that SNP are not linked
#'  make all GENETIC_POSITION = 0 (numeric).} 
#'  \item{'B_ALLELE_FREQ' is the frequency of the B allele in the population.} 
#'  \item{'ERROR_RATE' is the SNP error rate (i.e. the proportion of individuals/pools with signal intensity data 
#'  from a random genotype rather than the true genotype for the SNP_ID).  Refer to Hamilton 2020 (numeric).} 
#'  \item{'PROP_MISS' is the proportion missing data for the SNP_ID (numeric).} 
#' }

#' @param missing.parents is a vector idenifying parents with no SNP data (i.e. known missing parents). Samples/individuals in missing.parents must be present as a SIRE_ID or a DAM_ID in ped

#' @param true.snp.param.indiv is a data frame detailing the assumed SNP parameters of the population with the following 
#' headings (class in parentheses):
#' \itemize{
#'  \item{'SNP_ID' is the SNP identifier (character).} 
#'  \item{'MEAN_P_AA' is the mean of allelic proportion for homozygous A genotypes (numeric).}
#'  \item{'SD_P_AA' is the standard deviation of allelic proportion for homozygous A genotypes (numeric).} 
#'  \item{'MEAN_P_AB' is the mean of allelic proportion for heterozygous (AB) genotypes (numeric).}
#'  \item{'SD_P_AB' is the standard deviation of allelic proportion for heterozygous (AB) genotypes (numeric).}  
#'  \item{'MEAN_P_BB' is the mean of allelic proportion for homozygous B genotypes (numeric).}
#'  \item{'SD_P_BB' is the standard deviation of allelic proportion for homozygous B genotypes (numeric).} 
#'  \item{'A_ALLELE' is the base represented by allele A (i.e. 'A', 'C', 'G' or 'T') (character).}
#'  \item{'B_ALLELE' is the base represented by allele B (i.e. 'A', 'C', 'G' or 'T') (character).}
#' } 
#' @param sim.fam.sets is a data frame with the following headings (class in parentheses).  Note: if sim.fam.sets = NULL
#' (see example below with n.in.pools = 8), 
#' FAMILY_ID is taken from the 'fams' and duplicated n.in.pools times, FAM_SET_ID = 1 for the first duplication of 
#' FAMILY_IDs, 2 for the second etc and PROBABILITY = NA (default = NULL):
#' \itemize{
#'  \item{'FAM_SET_ID' is the family set identifier (integer).  A 'family set' is a group of families of which one is known to be the true family 
#'  of one of the individuals in a pooled sample.  Within each 'family set combination' there must be a 'family set'
#'  for each individual in a pooled sample (i.e. if n.in.pools = 2 there must be two family sets in each family set combination)} 
#'  \item{'FAMILY_ID' is the family identifier (integer).} 
#'  \item{'PROBABILITY' is probability that an individual from this family is represented in the pooled sample.  
#'  If all are NA it is assumed that the probability is equal for each family within the family set.} 
#' }
#' @param method is a vector of methods to be implemented (e.g. c("Quantitative", "Discrete", "Exclusion", "Least_squares"))
#' @param beta.min.ss is a logical variable appicable to least_squares method only (default = FALSE).
#' If TRUE, the sum of squares of all parental combinations are computed and the combination with the minimum value is identified.  
#' Refer to Hamilton 2020. 
#' @param discrete.method is a character variable applicable to the "Discrete" or "Exclusion" methods only
#' (default = "geno.probs").  It must equal either:
#' \itemize{
#'  \item{"geno.probs" in which case discrete genotypes for parents and pools are derived from genotype probabilities.}  
#'  \item{"assigned.genos" in which case discrete genotypes for parents and pools are obtained directly from the snp.dat.indiv and snp.dat.pools inputs.}
#'  }
#' @param threshold.indiv is a numeric variable between 0 and 1 inclusive applicable to the "Discrete" or "Exclusion" methods only
#' when discrete.method = "geno.probs" (default = NULL).  A discrete genotype is assigned to the the most likely genotype in 
#' the quantitative ordered genotype probability matrix Gij if it is greater than threshold.indiv (or
#' threshold.indiv / 2 for the two heterozygous genotypes).  Otherwise the genotype is deemed missing (refer to the left hand side of 
#' page 5 of Henshall et al. 2014)
#' @param threshold.pools is a numeric variable between 0 and 1 inclusive applicable to the "Discrete" or "Exclusion" methods only
#' when discrete.method = "geno.probs" (default = NULL).  Equivalent to threshold.indiv for pooled DNA samples.
#' @param n.in.pools is an integer variable representing the number of individual that contributed DNA to each sample in snp.dat.pools  
#' @param min.intensity is a numeric variable (default = 0).  If the square root of the sum of INTENSITY_A squared and 
#' INTENSITY_B squared in snp.dat.indiv or snp.dat.pools is less than min.intensity then this record is excluded.
#' That is, observations that fall into an arc with a radius equal to min.intensity in the lower left of
#' signal intensity scatter plots are excluded.  
#' @param snp.error.assumed Must be one of (default = NULL):
#' \itemize{
#'  \item{NULL.  Note that if snp.error.assumed is NULL then snp.error.underlying must not be NULL.}
#'  \item{a numeric variable between 0 and 1, in which case the 'assumed error rate' (see Henshall et al 2014) is the same across all SNP.}
#'  \item{a data frame with columns SNP_ID and SNP_ERROR_TILDE (see Henshall et al 2014).}
#' } 
#' @param snp.error.underlying. Not used if snp.error.assumed is not NULL (default = NULL). Must be either:
#' \itemize{
#'  \item{NULL.}
#'  \item{a numeric variable between 0 and 1 inclusive.  Used to comptute SNP_ERROR_TILDE from SNP_ERROR_HAT according
#'                      to the approach outlined on the left of page 5 of Henshall et al. 2014 using individual 
#'                      (i.e. not pooled) data only.  If snp.error.underlying = 0 then SNP_ERROR_TILDE = SNP_ERROR_HAT.}
#' } 
#' @param min.sd: a numeric variable defining a lower bound to be applied to estimates of the 
#' standard deviation of allelic proportion for genotypes in snp.param.indiv and snp.param.pools (default = 0)
#' @param fams is a data frame with the following headings (class in parentheses):
#' \itemize{
#'  \item{'FAMILY_ID' is the family identifier (integer).} 
#'  \item{'SIRE_ID' is the sire identifier (integer).} 
#'  \item{'DAM_ID' is the dam identifier (integer).} 
#' } 
#' @param skip.checks is a logical variable.  If FALSE parent.assign.fun data checks are not undertaken.


#' @return 'summary' is a data frame containing a summary of simulated pedigree assignments: 
#' \itemize{
#'  \item{'METHOD' is the method implemented.}           
#'  \item{'PARENTS_TO_ASSIGN' is a count of uncertain parents for which assignments were attempted.}
#'  \item{PROP_CORRECT_ASSIGN' is the proportion of PARENTS_TO_ASSIGN assigned correctly.}
#'  \item{'CRIT_DELTA_0.950' the delta LOD above which 95 percent of assignments were correct.  Applicable to 
#'  maximum likelihood methods only.}
#'  \item{'CRIT_DELTA_0.990' the delta LOD above which 99 percent of assignments were correct.  Applicable to 
#'  maximum likelihood methods only.}
#'  \item{'CRIT_DELTA_0.995' the delta LOD above which 99.5 percent of assignments were correct.  Applicable to 
#'  maximum likelihood methods only.}
#' } 
#' @return ggplot.log.quant: 
#' \itemize{
#'  \item{is a ggplot object (histogram) of delta LOD values using the 'Quantitative' method (if applicable).}
#' }
#' @return ggplot.log.discrete:
#' \itemize{
#'  \item{is a ggplot object (histogram) of delta LOD values using the Discrete' method (if applicable).} 
#' }
#' @return 'quant.sim.out' is a detailed summary for the 'Quantitative' method  (refer to Hamilton 2020):
#' \itemize{    
#'  \item{'TRUE_ID' is the true parent identifier}          
#'  \item{'REP' is the simualtion repetition number}
#'  \item{'PARENT_NUMBER' is a unique parent identifier within REP}
#'  \item{'QUANT_ID' is the assigned parent identifier using the 'Quantitative' method}
#'  \item{'QUANT_DELTA_LOD' is the delta LOD (refer to Hamilton 2020)}
#'  \item{'QUANT_LOD' is the LOD (refer to Hamilton 2020)}
#'  \item{'CORRECT_ASSIGN' is TRUE if the parent was correctly assigned}
#'  \item{'CUM_INCORRECT_ASSIGN' is a cumulative count of incorrectly assigned parents}
#'  \item{'CUM_PROP_CORRECT_ASSIGN' is cumulative proportion of correctly assigned parents}
#' } 
#' @return 'discrete.sim.out' is a detailed summary for the 'Discrete' method  (refer to Hamilton 2020):
#' \itemize{    
#'  \item{'TRUE_ID' is the true parent identifier.}           
#'  \item{'REP' is the simualtion repetition number.}
#'  \item{'PARENT_NUMBER' is a unique parent identifier within REP.}
#'  \item{'DISCRETE_ID' is the assigned parent identifier using the 'Discrete' method.}
#'  \item{'DISCRETE_DELTA_LOD' is the delta LOD.}
#'  \item{'DISCRETE_LOD' is the LOD.}
#'  \item{'CORRECT_ASSIGN' is TRUE if the parent was correctly assigned.}
#'  \item{'CUM_INCORRECT_ASSIGN' is a cumulative count of incorrectly assigned parents.}
#'  \item{'CUM_PROP_CORRECT_ASSIGN' is cumulative proportion of correctly assigned parents.}
#' } 
#' @return 'exclusion.sim.out' is a detailed summary for the 'Exclusion' method  (refer to Hamilton 2020): 
#' \itemize{    
#'  \item{'TRUE_ID' is the true parent identifier.}           
#'  \item{'REP' is the simualtion repetition number.}
#'  \item{'PARENT_NUMBER' is a unique parent identifier within REP.}
#'  \item{'EXCLUSION_ID' is the assigned parent identifier using the 'Exclusion' method.}
#'  \item{'CORRECT_ASSIGN' is TRUE if the parent was correctly assigned.}
#' } 
#' @return 'ls.sim.beta.constrain.out' is a detailed summary for the 'least squares' method where
#' beta hat is constrained to equal 1/n.in.pools within each FAM_SET_ID  (refer to Hamilton 2020): 
#' \itemize{    
#'  \item{'TRUE_ID' is the true parent identifier.}           
#'  \item{'REP' is the simualtion repetition number.}
#'  \item{'PARENT_NUMBER' is a unique parent identifier within REP.}
#'  \item{'LS_ID' is the assigned parent identifier using the 'least squares' method with beta hat constrained.}
#'  \item{'CORRECT_ASSIGN' is TRUE if the parent was correctly assigned.}
#' }  
#' @return 'ls.sim.min.ss.out' is a detailed summary for the 'least squares' method where
#' the family combination with the lowest sum of squares is identified  (refer to Hamilton 2020): 
#' \itemize{    
#'  \item{'TRUE_ID' is the true parent identifier.}          
#'  \item{'REP' is the simualtion repetition number.}
#'  \item{'PARENT_NUMBER' is a unique parent identifier within REP.}
#'  \item{'LS_ID' is the assigned parent identifier using the 'least squares' method with the lowest sum of squares is identified.}
#'  \item{'CORRECT_ASSIGN' is TRUE if the parent was correctly assigned.}
#' }  

#' @examples
#' #Retrieve data for 'pooling by phenotype' example from Hamilton 2020
#' data(shrimp.ped)
#' data(shrimp.map)
#' data(shrimp.true.snp.param.indiv)
#' data(shrimp.sim.fam.sets)
#' data(shrimp.fams)
#' 
#' #Run simulation for all methods with n.in.pools = 2.  Note that 3 is not enough repetitions (1000 may be).
#' sim.parent.assign.fun(n_repetitions = 3, 
#'                       ped = shrimp.ped,
#'                       map = shrimp.map,
#'                       true.snp.param.indiv = shrimp.true.snp.param.indiv,
#'                       sim.fam.sets = shrimp.sim.fam.sets, # equivalent to sim.fam.sets = NULL in this case
#'                       method = c("Quantitative", "Discrete", "Exclusion", "Least_squares"),     
#'                       beta.min.ss = TRUE, 
#'                       discrete.method = "geno.probs",   
#'                       threshold.indiv = 0.98,              
#'                       threshold.pools = 0.98, 
#'                       n.in.pools = 2,                
#'                       snp.error.assumed = 0.01,        
#'                       fams = shrimp.fams
#' )
#' 
#' #Run simulation using "Least_squares" method (beta.min.ss = FALSE) with n.in.pools = 8.  
#' #Do not attempt large pool sizes using any other method nor with beta.min.ss = TRUE, as your
#' #computer is likely to say no.  Note that 3 is not enough repetitions but is okay as an example.
#' sim.parent.assign.fun(n_repetitions = 3, 
#'                       ped = shrimp.ped,
#'                       map = shrimp.map,
#'                       true.snp.param.indiv = shrimp.true.snp.param.indiv,
#'                       sim.fam.sets = NULL, #shrimp.sim.fam.sets only appropriate for n.in.pools = 2
#'                       method = "Least_squares",     
#'                       beta.min.ss = FALSE,  
#'                       n.in.pools = 8,                
#'                       snp.error.assumed = 0.01,        
#'                       fams = shrimp.fams
#' )

#' @references Henshall JM, Dierens, L Sellars MJ (2014) Quantitative analysis of low-density SNP data for parentage assignment and estimation of family contributions to pooled samples. Genetics Selection Evolution 46, 51. https://doi 10.1186/s12711-014-0051-y 
#' @references Hamilton MG (2020) Maximum likelihood parentage assignment using quantitative genotypes


#' @export
 
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

sim.parent.assign.fun <- function(n_repetitions,
                                  ped,
                                  map,
                                  missing.parents = NULL, #vector of parents with no SNP data
                                  true.snp.param.indiv,
                                  sim.fam.sets = NULL,
                                  
                                  method, #c("Quantitative", "Discrete", "Exclusion", "Least_squares"),
                                  beta.min.ss = FALSE, #Appicable to least_squares method.  If TRUE, beta constrained to integers according to n.in.pools and minimum sum of squares identified
                                  discrete.method = "geno.probs", #"geno.probs" or "assigned.genos" 
                                  threshold.indiv = NULL,         #Dij.from.Gij.fun.  Not used if discrete.method = "assigned.genos"
                                  threshold.pools = NULL,         #dkj.from.gkj.fun.  Not used if discrete.method = "assigned.genos"
                                  
                                  #SNP data
                                  n.in.pools,
                                  min.intensity        = 0,    #pij.fun.  
                                  snp.error.assumed    = NULL, #If not null then this error is applied to all SNP.
                                  snp.error.underlying = NULL, #adj.geno.prob.fun.  Not required if snp.error.assumed is not NULL.
                                  
                                  #SNP parameters
                                  min.sd          = 0,         #phi.ij.fun   lambda.ij.fun  
                                  
                                  #Define pools
                                  fams,
                                  
                                  skip.checks = FALSE
) {
  print(Sys.time())
  print("Running sim.parent.assign.fun")  
  
  # load required packages
  if("dplyr" %in% installed.packages()[, "Package"] == FALSE) {install.packages("dplyr", repos='https://cran.csiro.au/')} 
  library(dplyr)
  if("reshape2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("reshape2", repos='https://cran.csiro.au/')} 
  library(reshape2)

  #Data checks####################################
  
  #n_repetitions checks
  n_repetitions <- as.integer(n_repetitions)
  if(n_repetitions < 1 ) {
    stop("n_repetitions must be an integer greater than 0")
  }
  
  #sim.fam.sets checks
  
  if(!is.null(sim.fam.sets)){
  
  #Check that required headings are present in sim.fam.sets  
  if(sum(c("FAM_SET_ID", "FAMILY_ID", "PROBABILITY") %in% colnames(sim.fam.sets)) != 3) {
    stop("sim.fam.sets input must be a data frame containing the following headings: FAM_SET_ID, FAMILY_ID, PROBABILITY")
  }
  
  sim.fam.sets <- sim.fam.sets[,c("FAM_SET_ID", "FAMILY_ID", "PROBABILITY")]
  
  sim.fam.sets$FAM_SET_ID  <- as.integer(sim.fam.sets$FAM_SET_ID)  
  sim.fam.sets$FAMILY_ID   <- as.integer(sim.fam.sets$FAMILY_ID)  
  sim.fam.sets$PROBABILITY <- as.numeric(sim.fam.sets$PROBABILITY)  
  
  if(sum(!is.integer(sim.fam.sets$FAM_SET_ID)) > 0) {
    stop("FAM_SET_ID in sim.fam.sets must be an integer.  Also check for missing values.")
  }
  if(sum(!is.integer(sim.fam.sets$FAMILY_ID)) > 0) {
    stop("FAMILY_ID in sim.fam.sets must be an integer.  Also check for missing values.")
  }
  if(sum(!is.numeric(sim.fam.sets$PROBABILITY)) > 0) {
    stop("PROBABILITY in sim.fam.sets must be numeric and between 0 and 1 (or NA).")
  }
  
  if(sum(sim.fam.sets[,c("FAM_SET_ID", "FAMILY_ID")] < 1) > 0) {
    stop("FAM_SET_ID and FAMILY_ID in sim.fam.sets must be greater than 0")
  }
  
  if((sum(sim.fam.sets[,"PROBABILITY"] > 1 & !is.na(sim.fam.sets[,"PROBABILITY"])) > 0) |
     (sum(sim.fam.sets[,"PROBABILITY"] < 0 & !is.na(sim.fam.sets[,"PROBABILITY"])) > 0)) {
    stop("PROBABILITY in sim.fam.sets must be between 0 and 1")
  }
  
  if((length(unique(sim.fam.sets[,"FAM_SET_ID"]))) != n.in.pools) {
    stop("The number of FAM_SET_IDs in sim.fam.sets does not equal n.in.pools")
  }
  
  } 
  
  #if sim.fam.sets is NULL then generate from fams
  if(is.null(sim.fam.sets)){
    sim.fam.sets <- data.frame(FAM_SET_ID = rep(1:n.in.pools, each = nrow(fams)),
      FAMILY_ID = rep(fams$FAMILY_ID,n.in.pools),
      PROBABILITY = NA)
  }                          
  
  #check probabilites and replace NAs if appropriate
  tmp.fam.sets <- unique(sim.fam.sets$FAM_SET_ID)
  for(tmp.fam.set in tmp.fam.sets) {
    if(sum(is.na(sim.fam.sets[sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set,"PROBABILITY"])) == 
       sum(sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set)) {
      sim.fam.sets[sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set,"PROBABILITY"] <- 
        1/sum(sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set) #change probabilities from NA to equal values within FAM_SET_ID
    } else {
      if(sum(sim.fam.sets[sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set,"PROBABILITY"], na.rm = T) != 1) {
        stop("PROBABILY within each level of FAM_SET_ID in sim.fam.sets must sum to 1 (or all be NA)")
      } else {
        sim.fam.sets[sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set &
                       is.na(sim.fam.sets[,"PROBABILITY"]),"PROBABILITY"] <- 0  #Replace NA PROBABILITY from sim.fam.sets with 0
      }
    }
  }
  
  #ped data checks
  if(sum(c("SAMPLE_ID", "SIRE_ID", "DAM_ID", "SAMPLED") %in% colnames(ped)) != 4) {
    stop("ped input must be a data frame containing the following headings: SAMPLE_ID, SIRE_ID, DAM_ID, SAMPLED")
  }
  
  ped <- ped[,c("SAMPLE_ID", "SIRE_ID", "DAM_ID", "SAMPLED")]
  
  ped[is.na(ped[,"DAM_ID"]),"DAM_ID"] <- 0
  ped[is.na(ped[,"SIRE_ID"]),"SIRE_ID"] <- 0
  
  ped$SAMPLE_ID  <- as.integer(ped$SAMPLE_ID)  
  ped$SIRE_ID  <- as.integer(ped$SIRE_ID)    
  ped$DAM_ID <- as.integer(ped$DAM_ID) 
  ped$SAMPLED <- as.logical(ped$SAMPLED) 
  
  if(0 %in% ped$SAMPLE_ID) {
    stop("SAMPLE_ID in ped must not be 0")
  }
  
  if(sum(!ped[,"DAM_ID"] %in% c(0,ped[,"SAMPLE_ID"])) > 0 |
     sum(!ped[,"SIRE_ID"] %in% c(0,ped[,"SAMPLE_ID"])) > 0) {
    stop("SIRE_IDs and DAM_IDs must be represented in the SAMPLE_ID column of ped")
  }
  
  if(length(ped[,"SAMPLE_ID"]) != nrow(ped)) {
    stop("Each SAMPLE_ID in ped must be unique")
  }
  
  if(sum(!((ped[,"DAM_ID"] > 0 & ped[,"SIRE_ID"] > 0) |
           (ped[,"DAM_ID"] == 0 & ped[,"SIRE_ID"] == 0))) > 0){
    stop("ped must consist of parents (i.e. both parents = 0) or full-sib individuals (i.e. both parents represented in the SAMPLE_ID column")
  }
  
  if(sum(is.na(ped$SAMPLED)) > 0 |
     !is.logical(ped$SAMPLED)) {
    stop("SAMPLED field in ped must be either TRUE or FALSE")
  }
  
  #map checks
  
  if(sum(c("CHROMOSOME",  "SNP_ID", "GENETIC_POSITION", "B_ALLELE_FREQ", "ERROR_RATE", "PROP_MISS") %in% colnames(map)) != 6) {
    stop("map input must be a data frame containing the following headings: CHROMOSOME, SNP_ID, GENETIC_POSITION, B_ALLELE_FREQ, ERROR_RATE, PROP_MISS")
  }
  
  map <- map[,c("CHROMOSOME",  "SNP_ID", "GENETIC_POSITION", "B_ALLELE_FREQ", "ERROR_RATE", "PROP_MISS")]
  
  map$CHROMOSOME  <- as.integer(map$CHROMOSOME)  
  map$SNP_ID  <- as.character(map$SNP_ID)    
  map$GENETIC_POSITION <- as.numeric(map$GENETIC_POSITION) 
  map$B_ALLELE_FREQ <- as.numeric(map$B_ALLELE_FREQ) 
  map$ERROR_RATE <- as.numeric(map$ERROR_RATE) 
  map$PROP_MISS <- as.numeric(map$PROP_MISS) 
  
  if(sum(map$GENETIC_POSITION < 0) > 0) {
    stop("GENETIC_POSITION in map must not be negative")
  }
  
  if(sum((map$B_ALLELE_FREQ < 0 | map$B_ALLELE_FREQ > 1)) > 0) {
    stop("B_ALLELE_FREQ in map must be between 0 and 1")
  }
  
  if(sum((map$ERROR_RATE < 0 | map$ERROR_RATE > 1)) > 0) {
    stop("ERROR_RATE in map must be between 0 and 1")
  }
  
  if(sum((map$PROP_MISS < 0 | map$PROP_MISS > 1)) > 0) {
    stop("PROP_MISS in map must be between 0 and 1")
  }
  
  #missing.parents checks
  
  if(!is.null(missing.parents)) {
    missing.parents <- as.integer(missing.parents)
    
    if(sum(!missing.parents %in% c(ped$SIRE_ID, ped$DAM_ID)) > 0) {
      stop("missing.parents must be a SIRE_ID or DAM_ID in ped")
    }
  }
  
  #true.snp.param.indiv checks
  
  if(sum(c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB", "A_ALLELE", "B_ALLELE") %in% colnames(true.snp.param.indiv)) != 9) {
    stop("true.snp.param.indiv input must be a data frame containing the following headings: SNP_ID MEAN_P_AA SD_P_AA MEAN_P_AB SD_P_AB MEAN_P_BB SD_P_BB A_ALLELE B_ALLELE")
  }
  
  true.snp.param.indiv <- true.snp.param.indiv[,c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB", "A_ALLELE", "B_ALLELE")]
  
  true.snp.param.indiv$SNP_ID           <- as.character(true.snp.param.indiv$SNP_ID)    
  true.snp.param.indiv$MEAN_P_AA <- as.numeric(true.snp.param.indiv$MEAN_P_AA) 
  true.snp.param.indiv$SD_P_AA    <- as.numeric(true.snp.param.indiv$SD_P_AA) 
  true.snp.param.indiv$MEAN_P_AB       <- as.numeric(true.snp.param.indiv$MEAN_P_AB) 
  true.snp.param.indiv$SD_P_AB        <- as.numeric(true.snp.param.indiv$SD_P_AB) 
  true.snp.param.indiv$MEAN_P_BB       <- as.numeric(true.snp.param.indiv$MEAN_P_BB) 
  true.snp.param.indiv$SD_P_BB        <- as.numeric(true.snp.param.indiv$SD_P_BB) 
  true.snp.param.indiv$A_ALLELE           <- as.character(true.snp.param.indiv$A_ALLELE) 
  true.snp.param.indiv$B_ALLELE           <- as.character(true.snp.param.indiv$B_ALLELE) 
  
  if(sum(true.snp.param.indiv[,"MEAN_P_AA"] > 1 | true.snp.param.indiv[,"MEAN_P_AA"] < 0 |
         true.snp.param.indiv[,"MEAN_P_AB"] > 1 | true.snp.param.indiv[,"MEAN_P_AB"] < 0 |
         true.snp.param.indiv[,"MEAN_P_BB"] > 1 | true.snp.param.indiv[,"MEAN_P_BB"] < 0 ) > 0) {
    stop("MEAN_P in true.snp.param.indiv must not be less than 0 or more than 1")
  }
  
  if(sum(true.snp.param.indiv[,"SD_P_AA"] < 0 |
         true.snp.param.indiv[,"SD_P_AB"] < 0 |
         true.snp.param.indiv[,"SD_P_BB"] < 0 ) > 0) {
    stop("SD_P in true.snp.param.indiv must not be less than 0")
  }
  
  #Check that A_ALLELE and B_ALLELE data are one of A, C, G, or T
  if(sum(true.snp.param.indiv[,"A_ALLELE"] %in% c("A", "C", "G", "T") & 
         true.snp.param.indiv[,"B_ALLELE"] %in% c("A", "C", "G", "T")) !=  nrow(true.snp.param.indiv)) {
    stop("A_ALLELE and B_ALLELE data in true.snp.param.indiv must be \'A\', \'C\', \'G\', or \'T\'")
  }
  
  #Check that A_ALLELE is not the same as B_ALLELE
  if(sum(true.snp.param.indiv[,"A_ALLELE"] == true.snp.param.indiv[,"B_ALLELE"]) > 0) {
    stop("A_ALLELE and B_ALLELE must be different for each SNP in true.snp.param.indiv")
  } 
  
  if(sum(!as.character(unique(map[,"SNP_ID"])) %in% as.character(unique(true.snp.param.indiv[,"SNP_ID"]))) != 0) {
    stop("not all SNP_ID in map are represented in true.snp.param.indiv")
  }
  
  #fam.set.combns checks
 # if(length(unique(fam.set.combns[,"FAM_SET_COMBN_ID"])) > 1) {
#    stop("For simulation there must only be one FAM_SET_COMBN_ID in fam.set.combns (i.e. FAM_SET_COMBN_ID must be the same in all rows")
#  }
  
#  fam.set.combns[,"FAM_SET_COMBN_ID"] <- 1
  
  #End data checks

  #create fam.set.combns
  fam.set.combns <- sim.fam.sets[,c("FAM_SET_ID", "FAMILY_ID")]
  fam.set.combns$FAM_SET_COMBN_ID <- 1
    
  #Seed data frames
  quant.sim.out <- NULL
  discrete.sim.out  <- NULL
  exclusion.sim.out <- NULL
  ls.sim.beta.constrain.out <- NULL
  ls.sim.min.ss.out <- NULL 
  
  #identify parents that must contribute to pool according to fam.set.combns
  tmp <- merge(fam.set.combns, fams, by = "FAMILY_ID", all.x = TRUE)
  fam.sets <- unique(fam.set.combns[,"FAM_SET_ID"])
  remove.parent <- NULL
  for(fam.set in fam.sets) {
    tmp.sires <- unique(tmp[tmp[,"FAM_SET_ID"] == fam.set,"SIRE_ID"])
    if(length(tmp.sires) == 1) {
      remove.parent <- c(remove.parent,tmp.sires)
    }
    tmp.dams <- unique(tmp[tmp[,"FAM_SET_ID"] == fam.set,"DAM_ID"])
    if(length(tmp.dams) == 1) {
      remove.parent <- c(remove.parent,tmp.dams)
    }
  }
  rm(tmp, fam.sets, tmp.sires, tmp.dams)
  
  ped.orig <- ped
  
  for(rep in 1:n_repetitions) {
  
    ped <- ped.orig  
    #create indivs.in.sim.pools - "FAM_SET_ID", "SAMPLE_ID"
    indivs.in.sim.pools <- NULL
    tmp.fam.sets <- unique(sim.fam.sets$FAM_SET_ID)
    for(tmp.fam.set in tmp.fam.sets) {
      tmp.sim.fam.sets <- sim.fam.sets[sim.fam.sets[,"FAM_SET_ID"] == tmp.fam.set,]
      tmp.indivs.in.sim.pools <- as.integer(cut(runif(1, 0, 1), breaks = c(0,cumsum(tmp.sim.fam.sets[,"PROBABILITY"]))))
      tmp.indivs.in.sim.pools <- data.frame(FAM_SET_ID = tmp.fam.set,
                                            FAMILY_ID = tmp.sim.fam.sets[tmp.indivs.in.sim.pools,"FAMILY_ID"])
      indivs.in.sim.pools <- rbind(indivs.in.sim.pools,tmp.indivs.in.sim.pools)
      rm(tmp.indivs.in.sim.pools, tmp.sim.fam.sets)
    }
    
    indivs.in.sim.pools <- merge(indivs.in.sim.pools, fams, by = "FAMILY_ID", all.x = TRUE)
    
    indivs.in.sim.pools$SAMPLE_ID <- c((max(ped$SAMPLE_ID)+1):(max(ped$SAMPLE_ID)+
                                                                 nrow(indivs.in.sim.pools)))
    indivs.in.sim.pools$SAMPLED <- FALSE
    ped <- rbind(ped[,c("SAMPLE_ID", "SIRE_ID", "DAM_ID","SAMPLED")],
                 indivs.in.sim.pools[,c("SAMPLE_ID", "SIRE_ID", "DAM_ID","SAMPLED")])
    
    indivs.in.sim.pools <- indivs.in.sim.pools[,c("FAM_SET_ID", "SAMPLE_ID")]
    
    #Get true parents
    true.sim.out.tmp <- ped[ped[,"SAMPLE_ID"] %in% indivs.in.sim.pools[,"SAMPLE_ID"],]
    true.sim.out.tmp <-   as.vector(unlist(true.sim.out.tmp[,c("SIRE_ID", "DAM_ID")]))
    true.sim.out.tmp <- true.sim.out.tmp[order(true.sim.out.tmp)]
    true.sim.out.tmp <- data.frame(TRUE_ID = true.sim.out.tmp,
                                   UNIQUE_TRUE_ID = paste(true.sim.out.tmp, 
                                                          sequence(rle(true.sim.out.tmp)$lengths), sep="_")
    )
    
    true.sim.out.tmp$REP <- rep
    
    #Get simulated data
    sim.pool <-  sim.pool.fun(map = map, 
                              ped = ped[,c("SAMPLE_ID", "SIRE_ID","DAM_ID")], 
                              true.snp.param.indiv = true.snp.param.indiv, 
                              indivs.in.sim.pools = indivs.in.sim.pools,
                              missing.parents = missing.parents)
    geno.indiv           <- sim.pool$geno.indiv
    geno.pool            <- sim.pool$geno.pool
    true.indivs.in.pool  <- sim.pool$true.indivs.in.pool
    sim.snp.dat.indiv    <- sim.pool$sim.snp.dat.indiv
    sim.snp.dat.pools     <- sim.pool$sim.snp.dat.pools
    true.snp.param.pools <- sim.pool$true.snp.param.pools
    rm(sim.pool)
    
    sim.snp.dat.pools$GENOTYPE <- sim.snp.dat.pools$OBS_GENO
    
    sim.snp.dat.indiv <- sim.snp.dat.indiv[sim.snp.dat.indiv[,"SAMPLE_ID"] %in% 
                                             ped[ped[,"SAMPLED"] == TRUE,"SAMPLE_ID"],]
    #Generate SNP parameters from simulated data
    sim.snp.param.indiv  <- snp.gen.param.fun(snp.dat.indiv = 
       sim.snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", 
                           "A_ALLELE", "B_ALLELE", "GENOTYPE")]) 
    
    sim.snp.param.pools <- snp.param.pools.fun(snp.param.indiv = sim.snp.param.indiv,
                                               n.in.pools=n.in.pools)
    
    fam.set.combns.by.pool <- data.frame(SAMPLE_ID = unique(geno.pool$SAMPLE_ID),
                                         FAM_SET_COMBN_ID = 1)
    
    #Assign parentage for simulated data
   
    running.sim <<- TRUE
    assignment <- parent.assign.fun(skip.checks = skip.checks,
                                    method = method,
                                    
                                    snp.dat.indiv = sim.snp.dat.indiv, 
                                    snp.dat.pools = sim.snp.dat.pools,
                                    min.intensity = min.intensity,             #pij.fun.  
                                    snp.error.assumed = snp.error.assumed,          #If not null then this error is applied to all SNP.
                                    snp.error.underlying = snp.error.underlying,          #adj.geno.prob.fun.  Not required if snp.error.assumed is not NULL.
                                    
                                    snp.param.indiv = sim.snp.param.indiv,
                                    snp.param.pools = sim.snp.param.pools,
                                    min.sd =min.sd,             #phi.ij.fun   lambda.ij.fun  
                                    
                                    fams = fams,
                                    n.in.pools = n.in.pools,
                                    fam.set.combns = fam.set.combns,
                                    fam.set.combns.by.pool = fam.set.combns.by.pool,
                                    
                                    beta.min.ss = beta.min.ss, #Appicable to least_squares method.  If TRUE, beta constrained to integers according to n.in.pools and minimum sum of squares identified
                                    
                                    discrete.method = discrete.method,      #"geno.probs" or "assigned.genos" 
                                    threshold.indiv = threshold.indiv,               #Dij.from.Gij.fun.  Not used if discrete.method = "assigned.genos"
                                    threshold.pools  = threshold.pools              #dkj.from.gkj.fun.  Not used if discrete.method = "assigned.genos"
                                    
    )
 
    rm(running.sim, pos = ".GlobalEnv")
    
    #get assigned parents and delta LOD 
    quant.sim.out.tmp      <- NULL
    discrete.sim.out.tmp   <- NULL
    exclusion.sim.out.tmp  <- NULL
    
    quant.delta.lod    <- NULL
    discrete.delta.lod <- NULL
    
    quant.lod    <- NULL
    discrete.lod <- NULL
    
    
    if("Quantitative" %in% method) {
      
      for (i in 1:(length(true.indivs.in.pool)*2)) {
        quant.sim.out.tmp      <- c(quant.sim.out.tmp, assignment$most.like.parents.quant[1,colnames(assignment$most.like.parents.quant)==paste("PARENT_", i, sep = "")])
        quant.delta.lod     <- c(quant.delta.lod, assignment$most.like.parents.quant[1,colnames(assignment$most.like.parents.quant)==paste("PARENT_", i, "_DELTA_LOD", sep = "")])
        quant.lod     <- c(quant.lod, assignment$most.like.parents.quant[1,colnames(assignment$most.like.parents.quant)=="LOD"])
      }
      
      quant.sim.out.tmp <- data.frame(PARENT_NUMBER = 1:length(quant.sim.out.tmp),
                                      QUANT_ID = quant.sim.out.tmp,
                                      UNIQUE_QUANT_ID = paste(quant.sim.out.tmp, 
                                                              sequence(rle(quant.sim.out.tmp)$lengths), sep="_"),
                                      QUANT_DELTA_LOD = quant.delta.lod,
                                      QUANT_LOD = quant.lod
      )
      quant.sim.out.tmp$CORRECT_ASSIGN <- quant.sim.out.tmp[,"UNIQUE_QUANT_ID"] %in% true.sim.out.tmp[,"UNIQUE_TRUE_ID"]
      quant.sim.out.tmp      <- quant.sim.out.tmp[order(quant.sim.out.tmp[,"CORRECT_ASSIGN"], decreasing = TRUE),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"]),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"] %in% quant.sim.out.tmp[,"UNIQUE_QUANT_ID"], decreasing = TRUE),]
      quant.sim.out.tmp     <- cbind(true.sim.out.tmp, quant.sim.out.tmp)
      #remove data for parents that must contribute to pool according to fam.set.combns
      quant.sim.out.tmp <- quant.sim.out.tmp[!quant.sim.out.tmp[,"TRUE_ID"] %in% remove.parent,]
      quant.sim.out         <- rbind(quant.sim.out, quant.sim.out.tmp)
    }
    
    if("Discrete" %in% method) {
      for (i in 1:(length(true.indivs.in.pool)*2)) {
        discrete.sim.out.tmp   <- c(discrete.sim.out.tmp, assignment$most.like.parents.discrete[1,colnames(assignment$most.like.parents.discrete)==paste("PARENT_", i, sep = "")])
        discrete.delta.lod  <- c(discrete.delta.lod, assignment$most.like.parents.discrete[1,colnames(assignment$most.like.parents.discrete)==paste("PARENT_", i, "_DELTA_LOD", sep = "")])
        discrete.lod  <- c(discrete.lod, assignment$most.like.parents.discrete[1,colnames(assignment$most.like.parents.discrete)=="LOD"])
      }
      
      discrete.sim.out.tmp  <- data.frame(PARENT_NUMBER = 1:length(discrete.sim.out.tmp),
                                          DISCRETE_ID = discrete.sim.out.tmp,
                                          UNIQUE_DISCRETE_ID = paste(discrete.sim.out.tmp, 
                                                                     sequence(rle(discrete.sim.out.tmp)$lengths), sep="_"),
                                          DISCRETE_DELTA_LOD = discrete.delta.lod,
                                          DISCRETE_LOD = discrete.lod
      )
      discrete.sim.out.tmp$CORRECT_ASSIGN <- discrete.sim.out.tmp[,"UNIQUE_DISCRETE_ID"] %in% true.sim.out.tmp[,"UNIQUE_TRUE_ID"]
      discrete.sim.out.tmp      <- discrete.sim.out.tmp[order(discrete.sim.out.tmp[,"CORRECT_ASSIGN"], decreasing = TRUE),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"]),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"] %in% discrete.sim.out.tmp[,"UNIQUE_DISCRETE_ID"], decreasing = TRUE),]
      discrete.sim.out.tmp  <- cbind(true.sim.out.tmp, discrete.sim.out.tmp)
      #remove data for parents that must contribute to pool according to fam.set.combns
      discrete.sim.out.tmp <- discrete.sim.out.tmp[!discrete.sim.out.tmp[,"TRUE_ID"] %in% remove.parent,]
      discrete.sim.out      <- rbind(discrete.sim.out, discrete.sim.out.tmp)
    }
    
    if("Exclusion" %in% method) {
      for (i in 1:(length(true.indivs.in.pool)*2)) {
        exclusion.sim.out.tmp  <- c(exclusion.sim.out.tmp, assignment$most.like.parents.excl.non.dup[1,colnames(assignment$most.like.parents.excl.non.dup)==paste("PARENT_", i, sep = "")])
      }
      
      exclusion.sim.out.tmp <- data.frame(PARENT_NUMBER = 1:length(exclusion.sim.out.tmp),
                                          EXCLUSION_ID = exclusion.sim.out.tmp,
                                          UNIQUE_EXCLUSION_ID = paste(exclusion.sim.out.tmp, 
                                                                      sequence(rle(exclusion.sim.out.tmp)$lengths), sep="_")
      )
      exclusion.sim.out.tmp$CORRECT_ASSIGN <- exclusion.sim.out.tmp[,"UNIQUE_EXCLUSION_ID"] %in% true.sim.out.tmp[,"UNIQUE_TRUE_ID"]
      exclusion.sim.out.tmp      <- exclusion.sim.out.tmp[order(exclusion.sim.out.tmp[,"CORRECT_ASSIGN"], decreasing = TRUE),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"]),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"] %in% exclusion.sim.out.tmp[,"UNIQUE_EXCLUSION_ID"], decreasing = TRUE),]
      exclusion.sim.out.tmp <- cbind(true.sim.out.tmp, exclusion.sim.out.tmp)
      #remove data for parents that must contribute to pool according to fam.set.combns
      exclusion.sim.out.tmp <- exclusion.sim.out.tmp[!exclusion.sim.out.tmp[,"TRUE_ID"] %in% remove.parent,]
      exclusion.sim.out     <- rbind(exclusion.sim.out, exclusion.sim.out.tmp)
    }
    
    #ls with beta.min.ss
    
    if("Least_squares" %in% method) {
      
      tmp <- assignment$beta[assignment$beta[,"BETA_HAT_CONSTRAINED"] != 0,]
      tmp$BETA_HAT_CONSTRAINED <- as.integer(tmp$BETA_HAT_CONSTRAINED * n.in.pools)
      ls.sim.beta.constrain.out.tmp  <- c(rep(tmp$SIRE_ID,tmp$BETA_HAT_CONSTRAINED), 
                                          rep(tmp$DAM_ID,tmp$BETA_HAT_CONSTRAINED)) #list of parent (duplicates included)
      ls.sim.beta.constrain.out.tmp  <- as.numeric(ls.sim.beta.constrain.out.tmp[order(ls.sim.beta.constrain.out.tmp)])
      rm(tmp)

       ls.sim.beta.constrain.out.tmp  <- data.frame(PARENT_NUMBER = 1:length(ls.sim.beta.constrain.out.tmp),
                                                   LS_ID = ls.sim.beta.constrain.out.tmp,
                                                   UNIQUE_LS_ID = paste(ls.sim.beta.constrain.out.tmp, 
                                                                        sequence(rle(ls.sim.beta.constrain.out.tmp)$lengths), sep="_")
      )
      
      ls.sim.beta.constrain.out.tmp$CORRECT_ASSIGN <- ls.sim.beta.constrain.out.tmp[,"UNIQUE_LS_ID"] %in% true.sim.out.tmp[,"UNIQUE_TRUE_ID"]
      ls.sim.beta.constrain.out.tmp      <- ls.sim.beta.constrain.out.tmp[order(ls.sim.beta.constrain.out.tmp[,"CORRECT_ASSIGN"], decreasing = TRUE),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"]),]
      true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"] %in% ls.sim.beta.constrain.out.tmp[,"UNIQUE_LS_ID"], decreasing = TRUE),]
      ls.sim.beta.constrain.out.tmp    <- cbind(true.sim.out.tmp, ls.sim.beta.constrain.out.tmp)
      #remove data for parents that must contribute to pool according to fam.set.combns
      ls.sim.beta.constrain.out.tmp <- ls.sim.beta.constrain.out.tmp[!ls.sim.beta.constrain.out.tmp[,"TRUE_ID"] %in% remove.parent,]
      ls.sim.beta.constrain.out        <- rbind(ls.sim.beta.constrain.out, ls.sim.beta.constrain.out.tmp)
      
      
      if(beta.min.ss) {
        tmp <- assignment$beta[assignment$beta[,"BETA_MIN_SS"] != 0,]
        
        ls.sim.min.ss.out.tmp  <- c(tmp$SIRE_ID, tmp$DAM_ID)
        ls.sim.min.ss.out.tmp  <- as.numeric(ls.sim.min.ss.out.tmp[order(ls.sim.min.ss.out.tmp)])
        rm(tmp)
  
        ls.sim.min.ss.out.tmp  <- data.frame(PARENT_NUMBER = 1:length(ls.sim.min.ss.out.tmp),
                                             LS_ID = ls.sim.min.ss.out.tmp,
                                             UNIQUE_LS_ID = paste(ls.sim.min.ss.out.tmp, 
                                                                  sequence(rle(ls.sim.min.ss.out.tmp)$lengths), sep="_")
        )
        
        ls.sim.min.ss.out.tmp$CORRECT_ASSIGN <- ls.sim.min.ss.out.tmp[,"UNIQUE_LS_ID"] %in% true.sim.out.tmp[,"UNIQUE_TRUE_ID"]
        ls.sim.min.ss.out.tmp      <- ls.sim.min.ss.out.tmp[order(ls.sim.min.ss.out.tmp[,"CORRECT_ASSIGN"], decreasing = TRUE),]
        true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"]),]
        true.sim.out.tmp      <- true.sim.out.tmp[order(true.sim.out.tmp[,"UNIQUE_TRUE_ID"] %in% ls.sim.min.ss.out.tmp[,"UNIQUE_LS_ID"], decreasing = TRUE),]
        ls.sim.min.ss.out.tmp    <- cbind(true.sim.out.tmp, ls.sim.min.ss.out.tmp)
        #remove data for parents that must contribute to pool according to fam.set.combns
        ls.sim.min.ss.out.tmp <- ls.sim.min.ss.out.tmp[!ls.sim.min.ss.out.tmp[,"TRUE_ID"] %in% remove.parent,]
        ls.sim.min.ss.out        <- rbind(ls.sim.min.ss.out, ls.sim.min.ss.out.tmp)
      } else {
        ls.sim.min.ss.out <- data.frame(CORRECT_ASSIGN = NA)
      }
    }
    
    print("###############################################################################################")
    print(paste("End Rep", rep, "of", n_repetitions, Sys.time()))
    print("###############################################################################################")
    
  } #end for(rep in 1:n_repetitions) {

  #Summarise and plot outputs
  
  n.quant.sim <- nrow(quant.sim.out)
  n.discrete.sim <- nrow(discrete.sim.out)
  n.exclusion.sim <- nrow(exclusion.sim.out)
  n.ls.sim.beta.constrain <- nrow(ls.sim.beta.constrain.out)
  n.ls.sim.min.ss <- nrow(ls.sim.min.ss.out) 
  
  if(is.null(n.quant.sim))  {n.quant.sim <- 0}
  if(is.null(n.discrete.sim))  {n.discrete.sim <- 0}
  if(is.null(n.exclusion.sim))  {n.exclusion.sim <- 0}
  if(is.null(n.ls.sim.beta.constrain))  {n.ls.sim.beta.constrain <- 0}
  if(is.null(n.ls.sim.min.ss))  {n.ls.sim.min.ss <- 0}
  if(ncol(ls.sim.min.ss.out) == 1)  {n.ls.sim.min.ss <- 0}
  
  prop.corr.quant.sim <- sum(quant.sim.out[,"CORRECT_ASSIGN"]) / n.quant.sim #quantitative proportion TRUE
  prop.corr.discrete.sim <- sum(discrete.sim.out[,"CORRECT_ASSIGN"]) / n.discrete.sim #discrete proportion TRUE
  prop.corr.exclusion.sim <-  sum(exclusion.sim.out[,"CORRECT_ASSIGN"]) / n.exclusion.sim #exclusion proportion TRUE
  prop.corr.ls.sim.beta.constrain <- sum(ls.sim.beta.constrain.out[,"CORRECT_ASSIGN"]) / n.ls.sim.beta.constrain
  prop.corr.ls.sim.min.ss <-  sum(ls.sim.min.ss.out[,"CORRECT_ASSIGN"]) / n.ls.sim.min.ss 
  
  #Get quantitative parameters
  quant.crit.delta.950 <- NA
  quant.crit.delta.990 <- NA
  quant.crit.delta.995 <- NA
  
  if("Quantitative" %in% method) {
    
    #get Critical Delta
    quant.sim.out <- quant.sim.out[order(quant.sim.out$QUANT_DELTA_LOD, decreasing = TRUE),]
    tmp_correct <- cumsum(quant.sim.out$CORRECT_ASSIGN) 
    quant.sim.out$CUM_INCORRECT_ASSIGN <- cumsum(!quant.sim.out$CORRECT_ASSIGN) 
    quant.sim.out$CUM_PROP_CORRECT_ASSIGN <- tmp_correct / (tmp_correct + quant.sim.out$CUM_INCORRECT_ASSIGN)
    rm(tmp_correct)
    
    if(min(quant.sim.out$CUM_PROP_CORRECT_ASSIGN) < 0.995) {
      tmp <- min(quant.sim.out[quant.sim.out$CUM_PROP_CORRECT_ASSIGN < 0.995,"CUM_INCORRECT_ASSIGN"])-1
      quant.crit.delta.995 <- quant.sim.out[quant.sim.out[,"CUM_INCORRECT_ASSIGN"] == tmp & !quant.sim.out[,"CORRECT_ASSIGN"],"QUANT_DELTA_LOD"]
      if(length(quant.crit.delta.995) == 0) {quant.crit.delta.995 <- NA}
      rm(tmp)
    }
    
    if(min(quant.sim.out$CUM_PROP_CORRECT_ASSIGN) < 0.990) {
      tmp <- min(quant.sim.out[quant.sim.out$CUM_PROP_CORRECT_ASSIGN < 0.990,"CUM_INCORRECT_ASSIGN"])-1
      quant.crit.delta.990 <- quant.sim.out[quant.sim.out[,"CUM_INCORRECT_ASSIGN"] == tmp & !quant.sim.out[,"CORRECT_ASSIGN"],"QUANT_DELTA_LOD"]
      if(length(quant.crit.delta.990) == 0) {quant.crit.delta.990 <- NA}
      rm(tmp)
    }
    
    if(min(quant.sim.out$CUM_PROP_CORRECT_ASSIGN) < 0.950) {
      tmp <- min(quant.sim.out[quant.sim.out$CUM_PROP_CORRECT_ASSIGN < 0.950,"CUM_INCORRECT_ASSIGN"])-1
      quant.crit.delta.950 <- quant.sim.out[quant.sim.out[,"CUM_INCORRECT_ASSIGN"] == tmp & !quant.sim.out[,"CORRECT_ASSIGN"],"QUANT_DELTA_LOD"]
      if(length(quant.crit.delta.950) == 0) {quant.crit.delta.950 <- NA}
      rm(tmp)
    }
    
  }
  
  #Get discrete parameters
  discrete.crit.delta.950 <- NA
  discrete.crit.delta.990 <- NA   
  discrete.crit.delta.995 <- NA 
  
  if("Discrete" %in% method) {

    #get Critical Delta
    discrete.sim.out <- discrete.sim.out[order(discrete.sim.out$DISCRETE_DELTA_LOD, decreasing = TRUE),]
    tmp_correct <- cumsum(discrete.sim.out$CORRECT_ASSIGN) 
    discrete.sim.out$CUM_INCORRECT_ASSIGN <- cumsum(!discrete.sim.out$CORRECT_ASSIGN) 
    discrete.sim.out$CUM_PROP_CORRECT_ASSIGN <- tmp_correct / (tmp_correct + discrete.sim.out$CUM_INCORRECT_ASSIGN)
    rm(tmp_correct)
    
    if(min(discrete.sim.out$CUM_PROP_CORRECT_ASSIGN) < 0.995) {
      tmp <- min(discrete.sim.out[discrete.sim.out$CUM_PROP_CORRECT_ASSIGN < 0.995,"CUM_INCORRECT_ASSIGN"])-1
      discrete.crit.delta.995 <- discrete.sim.out[discrete.sim.out[,"CUM_INCORRECT_ASSIGN"] == tmp & !discrete.sim.out[,"CORRECT_ASSIGN"],"DISCRETE_DELTA_LOD"]
      if(length(discrete.crit.delta.995) == 0) {discrete.crit.delta.995 <- NA}
      rm(tmp)
    }
    
    if(min(discrete.sim.out$CUM_PROP_CORRECT_ASSIGN) < 0.990) {
      tmp <- min(discrete.sim.out[discrete.sim.out$CUM_PROP_CORRECT_ASSIGN < 0.990,"CUM_INCORRECT_ASSIGN"])-1
      discrete.crit.delta.990 <- discrete.sim.out[discrete.sim.out[,"CUM_INCORRECT_ASSIGN"] == tmp & !discrete.sim.out[,"CORRECT_ASSIGN"],"DISCRETE_DELTA_LOD"]
      if(length(discrete.crit.delta.990) == 0) {discrete.crit.delta.990 <- NA}
      rm(tmp)
    }
    
    if(min(discrete.sim.out$CUM_PROP_CORRECT_ASSIGN) < 0.950) {
      tmp <- min(discrete.sim.out[discrete.sim.out$CUM_PROP_CORRECT_ASSIGN < 0.950,"CUM_INCORRECT_ASSIGN"])-1
      discrete.crit.delta.950 <- discrete.sim.out[discrete.sim.out[,"CUM_INCORRECT_ASSIGN"] == tmp & !discrete.sim.out[,"CORRECT_ASSIGN"],"DISCRETE_DELTA_LOD"]
      if(length(discrete.crit.delta.950) == 0) {discrete.crit.delta.950 <- NA}
      rm(tmp)
    }
  }
  
  summary <- data.frame(METHOD = c("Quantitative ML", "Discrete ML", "Exclusion", "Least squares beta hat constrained", "Least squares min SS"),
                        PARENTS_TO_ASSIGN = c(n.quant.sim,
                                              n.discrete.sim,
                                              n.exclusion.sim,
                                              n.ls.sim.beta.constrain,
                                              n.ls.sim.min.ss 
                        ),
                        PROP_CORRECT_ASSIGN = c(prop.corr.quant.sim,
                                                prop.corr.discrete.sim ,
                                                prop.corr.exclusion.sim, 
                                                prop.corr.ls.sim.beta.constrain ,
                                                prop.corr.ls.sim.min.ss  
                        ),
                        CRIT_DELTA_0.950 = c(quant.crit.delta.950, 
                                             discrete.crit.delta.950, #discrete maximum LOD of individual with FALSE assignment
                                             NA,
                                             NA,
                                             NA),
                        CRIT_DELTA_0.990 = c(quant.crit.delta.990, 
                                             discrete.crit.delta.990, #discrete maximum LOD of individual with FALSE assignment
                                             NA,
                                             NA,
                                             NA),
                        CRIT_DELTA_0.995 = c(quant.crit.delta.995, #quantitative maximum LOD of individual with FALSE assignment
                                             discrete.crit.delta.995, #discrete maximum LOD of individual with FALSE assignment
                                             NA,
                                             NA,
                                             NA)
                        
  )
  
  summary$METHOD <- as.character(summary$METHOD)
  summary <- summary[summary[,"PARENTS_TO_ASSIGN"] != 0,]
  
  #ls proportion TRUE
  
  library(ggplot2)
  
  bw <- 0.4
  
  #plots
  ggplot.log.quant <- NULL
  if("Quantitative" %in% method) {
    dat <- quant.sim.out[,c("CORRECT_ASSIGN","QUANT_DELTA_LOD")]
    dat$CORRECT_ASSIGN <- as.factor(dat$CORRECT_ASSIGN)
    #dat[is.na(dat[,"QUANT_DELTA_LOD"]),"QUANT_DELTA_LOD"] <- 0
    ggplot.log.quant <- ggplot(dat,aes(x=QUANT_DELTA_LOD)) + 
      geom_histogram(data=subset(dat,CORRECT_ASSIGN == 'TRUE'),fill = "grey10", alpha = 0.2, binwidth = bw) +#
      stat_bin(data=subset(dat,CORRECT_ASSIGN == 'FALSE'), geom = 'step', binwidth = bw, position=position_nudge(x=-bw/2))  +  #, alpha = 0.2
      ylab("Frequency") +
      xlab(expression(paste(Delta,"LOD"))) +
      theme_bw()  
  }
  
  ggplot.log.discrete <- NULL
  if("Discrete" %in% method) {
    dat <- discrete.sim.out[,c("CORRECT_ASSIGN","DISCRETE_DELTA_LOD")]
    dat$CORRECT_ASSIGN <- as.factor(dat$CORRECT_ASSIGN)
    #dat[is.na(dat[,"DISCRETE_DELTA_LOD"]),"DISCRETE_DELTA_LOD"] <- 0
    ggplot.log.discrete <- ggplot(dat,aes(x=DISCRETE_DELTA_LOD)) + 
      geom_histogram(data=subset(dat,CORRECT_ASSIGN == 'TRUE'),fill = "grey10", alpha = 0.2, binwidth = bw) +#
      stat_bin(data=subset(dat,CORRECT_ASSIGN == 'FALSE'), geom = 'step', binwidth = bw, position=position_nudge(x=-bw/2))  +  #, alpha = 0.2
      ylab("Frequency") +
      xlab(expression(paste(Delta,"LOD"))) + 
      theme_bw()  
  }
  
  #remove UNIQUE_TRUE_ID and UNIQUE_QUANT_ID and equivalents from quant.sim.out and discrete.sim.out
  quant.sim.out <-  quant.sim.out[,!colnames(quant.sim.out) %in% c("UNIQUE_TRUE_ID", "UNIQUE_QUANT_ID")]
  discrete.sim.out <-  discrete.sim.out[,!colnames(discrete.sim.out) %in% c("UNIQUE_TRUE_ID", "UNIQUE_DISCRETE_ID")]
  exclusion.sim.out <-  exclusion.sim.out[,!colnames(exclusion.sim.out) %in% c("UNIQUE_TRUE_ID", "UNIQUE_EXCLUSION_ID")]
  ls.sim.beta.constrain.out <-  ls.sim.beta.constrain.out[,!colnames(ls.sim.beta.constrain.out) %in% c("UNIQUE_TRUE_ID", "UNIQUE_LS_ID")]
  ls.sim.min.ss.out <-  ls.sim.min.ss.out[,!colnames(ls.sim.min.ss.out) %in% c("UNIQUE_TRUE_ID", "UNIQUE_LS_ID")]

  return(list(quant.sim.out     = quant.sim.out, 
              discrete.sim.out  = discrete.sim.out, 
              exclusion.sim.out = exclusion.sim.out, 
              ls.sim.beta.constrain.out = ls.sim.beta.constrain.out,
              ls.sim.min.ss.out = ls.sim.min.ss.out, 
              summary = summary,
              ggplot.log.quant = ggplot.log.quant,
              ggplot.log.discrete = ggplot.log.discrete))
  
} 

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

sim.geno.indiv.fun <- function(map, ped) {
  
  #map:
  #CHROMOSOME (integer), SNP_ID (character), GENETIC_POSITION (Morgans - numeric), 
  #  B_ALLELE_FREQ (numeric - 0-1), ERROR_RATE (numeric 0-1) , PROP_MISS (numeric 0-1)
  
  #(ordered by physical position within chromosome)
  #if wish to assume all SNP are independent (i.e. not linked) give all SNP different CHROMOSOME IDs and make all GENETIC_POSITION = 0
  
  #ped:
  #SAMPLE_ID (integer), SIRE_ID (integer), DAM_ID (integer)
  #SIRE_ID and DAM_ID should be 0 if unknown.  Parents with no SAMPLE_ID should be assigned a dummy SAMPLE_ID.
  
  map[,"SNP_ID"]   <- as.character(map[,"SNP_ID"])
  
  #Check column names
  
  #Get genetic distance within chromosomes
  map[,"GENETIC_DISTANCE"] <- c(0,diff(map[,"GENETIC_POSITION"]))
  map[c(0,diff(map[,"CHROMOSOME"])) == 1,"GENETIC_DISTANCE"] <- 0
  
  #define function to identify even numbers
  is.even.fun <- function(x) x %% 2 == 0
  
  #cycle through samples
  geno <- NULL
  for(i in 1:nrow(ped)) {
    # for(i in 1:length(founders)) {  
    sample <- ped[i,"SAMPLE_ID"]
    sire <- ped[ped[,"SAMPLE_ID"] == sample,"SIRE_ID"]
    dam  <- ped[ped[,"SAMPLE_ID"] == sample,"DAM_ID"]
    
    #get sire gamete haplo
    if(sire == 0) {
      sire.haplo <- map[,c("SNP_ID","B_ALLELE_FREQ")]
      sire.haplo[,"TRUE_SIRE_ALLELE"] <- "A"
      sire.haplo[runif(n = nrow(sire.haplo)) < sire.haplo[,"B_ALLELE_FREQ"],"TRUE_SIRE_ALLELE"] <- "B"
      #   sire.haplo <- sire.haplo[,colnames(sire.haplo) != "B_ALLELE_FREQ"]
    } else {
      sire.geno <- geno[geno[,"SAMPLE_ID"] == sire,]
      #identify cross over positions
      cross.over <- runif(n = nrow(map)) < map[,"GENETIC_DISTANCE"]
      #cross.over random for first position in each chromosome
      cross.over[c(1,diff(map[,"CHROMOSOME"])) == 1] <- runif(n = sum(c(1,diff(map[,"CHROMOSOME"])) == 1)) < 0.5
      sire.haplo <- sire.geno[,c("SNP_ID", "B_ALLELE_FREQ", "TRUE_SIRE_ALLELE")]
      use.dam.allele <- is.even.fun(cumsum(cross.over))
      sire.haplo[use.dam.allele,"TRUE_SIRE_ALLELE"] <- sire.geno[use.dam.allele,"TRUE_DAM_ALLELE"]
      rm(sire.geno, cross.over, use.dam.allele)
    }
    
    #get dam gamete haplo
    if(dam == 0) {
      dam.haplo <- map[,c("SNP_ID","B_ALLELE_FREQ")]
      dam.haplo[,"TRUE_DAM_ALLELE"] <- "A"
      dam.haplo[runif(n = nrow(dam.haplo)) < dam.haplo[,"B_ALLELE_FREQ"],"TRUE_DAM_ALLELE"] <- "B"
      #  dam.haplo <- dam.haplo[,colnames(dam.haplo) != "B_ALLELE_FREQ"]
    } else {
      dam.geno <- geno[geno[,"SAMPLE_ID"] == dam,]
      #identify cross over positions
      cross.over <- runif(n = nrow(map)) < map[,"GENETIC_DISTANCE"]
      #cross.over random for first position in each chromosome
      cross.over[c(1,diff(map[,"CHROMOSOME"])) == 1] <- runif(n = sum(c(1,diff(map[,"CHROMOSOME"])) == 1)) < 0.5
      dam.haplo <- dam.geno[,c("SNP_ID", "B_ALLELE_FREQ", "TRUE_DAM_ALLELE")]
      use.sire.allele <- is.even.fun(cumsum(cross.over))
      dam.haplo[use.sire.allele,"TRUE_DAM_ALLELE"] <- dam.geno[use.sire.allele,"TRUE_DAM_ALLELE"]
      rm(dam.geno, cross.over, use.sire.allele)
    }
    
    samp.geno <- merge(sire.haplo, dam.haplo, by = c("SNP_ID", "B_ALLELE_FREQ"), all = TRUE, sort = FALSE)
    
    #Use genotype replacement model of Marshall et al (1998)
    snp.id.error <-  map[runif(n = nrow(map)) < map[,"ERROR_RATE"],"SNP_ID"]  
    
    samp.geno[,"GENO_ERROR"] <- FALSE  
    samp.geno[,"OBS_SIRE_ALLELE"] <- samp.geno[,"TRUE_SIRE_ALLELE"]
    samp.geno[,"OBS_DAM_ALLELE"] <- samp.geno[,"TRUE_DAM_ALLELE"]
    
    if(length(snp.id.error) > 0) {
      snp.error <- samp.geno[,"SNP_ID"] %in% snp.id.error
      
      samp.geno[snp.error,"GENO_ERROR"] <- TRUE
      samp.geno[snp.error,"OBS_SIRE_ALLELE"] <- "A"
      replace <-  runif(n = nrow(samp.geno)) < samp.geno[,"B_ALLELE_FREQ"]
      samp.geno[snp.error &
                  replace ,
                "OBS_SIRE_ALLELE"] <- rep("B",sum(snp.error & replace))
      
      samp.geno[snp.error,"OBS_DAM_ALLELE"] <- "A"
      replace <-  runif(n = nrow(samp.geno)) < samp.geno[,"B_ALLELE_FREQ"]
      samp.geno[snp.error &
                  replace,
                "OBS_DAM_ALLELE"] <- rep("B",sum(snp.error & replace))
      rm(replace, snp.error)
    }
    
    samp.geno[,"SAMPLE_ID"] <- sample
    
    #retain only necessary columns
    samp.geno <- samp.geno[,c("SNP_ID", "B_ALLELE_FREQ", "TRUE_SIRE_ALLELE","TRUE_DAM_ALLELE", "SAMPLE_ID", "GENO_ERROR", "OBS_SIRE_ALLELE","OBS_DAM_ALLELE")]
    
    #rbind with previous SAMPLE_ID
    geno <- rbind(geno, samp.geno)
  }
  
  #add "A_ALLELE", "B_ALLELE"
  #  geno <- merge(geno, map[,c("SNP_ID", "A_ALLELE", "B_ALLELE")], by = "SNP_ID", all.x = TRUE)
  
  #get true concatenated unordered genotype
  geno[,"TRUE_GENO"] <- paste(geno[,"TRUE_SIRE_ALLELE"], geno[,"TRUE_DAM_ALLELE"], sep = "")
  geno[geno[,"TRUE_GENO"] == "BA","TRUE_GENO"] <- "AB" 
  
  #get observed concatenated unordered genotype
  geno[,"OBS_GENO"] <- paste(geno[,"OBS_SIRE_ALLELE"], geno[,"OBS_DAM_ALLELE"], sep = "")
  geno[geno[,"OBS_GENO"] == "BA","OBS_GENO"] <- "AB" 
  
  geno <- geno[,c("SNP_ID", "B_ALLELE_FREQ", "TRUE_GENO", "SAMPLE_ID", "GENO_ERROR", "OBS_GENO")]
  
  return(geno)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

sim.snp.dat.indiv.fun <- function(geno, true.snp.param.indiv) {
  
  #geno: output of sim.geno.fun
  
  # true.snp.param.indiv: Data frame. Output of snp.gen.param.fun
  #              SNP_ID        is the SNP identifier, 
  #              MEAN_P_AA     is the mean of allelic proportion (homozygous allele A), 
  #              SD_P_AA       is the standard deviation of allelic proportion (homozygous allele A), 
  #              MEAN_P_AB     is the mean of allelic proportion (heterozygous), 
  #              SD_P_AB       is the standard deviation of allelic proportion (heterozygous), 
  #              MEAN_P_BB     is the mean of allelic proportion (homozygous allele B), 
  #              SD_P_BB       is the standard deviation of allelic proportion (homozygous allele B),
  #              A_ALLELE is the base represented by the A allele (A, C, G, T)
  #              B_ALLELE is the base represented by the B allele (A, C, G, T)
  
  # snp.dat.indiv: Data frame (individuals not included present as a sire or dam in fams are considered offspring)
  #              SAMPLE_ID is the individual identifier
  #              SNP_ID   is the SNP identifier
  #              INTENSITY_A   is the area/intensity for allele A
  #              INTENSITY_B   is the area/intensity for allele B
  #              A_ALLELE  is the base represented by the A allele
  #              B_ALLELE  is the base represented by the B allele
  #              GENOTYPE is the SNP genotype call
  
  true.snp.param.indiv[,"SNP_ID"]   <- as.character(true.snp.param.indiv[,"SNP_ID"])
  true.snp.param.indiv[,"A_ALLELE"] <- as.character(true.snp.param.indiv[,"A_ALLELE"])
  true.snp.param.indiv[,"B_ALLELE"] <- as.character(true.snp.param.indiv[,"B_ALLELE"])
  
  geno[,"TRUE_GENO"] <- as.character(geno[,"TRUE_GENO"])
  geno[,"OBS_GENO"] <- as.character(geno[,"OBS_GENO"])
  
  n.in.pools <- nchar(geno[1,"TRUE_GENO"])/2
  
  #Get allelic proportion
  #geno <- merge(geno, true.snp.param.indiv, by = "SNP_ID", all.x = TRUE)
  geno$SNP_ID    <- as.character(geno$SNP_ID)
  true.snp.param.indiv$SNP_ID    <- as.character(true.snp.param.indiv$SNP_ID)
  geno <- left_join(geno, true.snp.param.indiv, by = "SNP_ID")
  
  geno <- geno[order(geno[,"SAMPLE_ID"]),]
  
  genotypes <- genotypes.fun(n.in.pools*2)
  
  for(genotype in genotypes) {
    tmp.obs <-  geno[,"OBS_GENO"] == genotype
    tmp.mean <- geno[,paste("MEAN_P_", genotype, sep = "")]
    tmp.sd   <- geno[,paste("SD_P_", genotype, sep = "")]
    geno[tmp.obs,"OBS_ALLELIC_PROP"] <- 
      rnorm(n = sum(tmp.obs)) * tmp.sd[tmp.obs] + tmp.mean[tmp.obs] 
    rm(tmp.obs, tmp.mean, tmp.sd)
  }
  
  #convert genotypes to bases
  for(snp in unique(geno$SNP_ID)) {
    a.allele <- geno[geno[,"SNP_ID"] == snp,"A_ALLELE"][1]
    b.allele <- geno[geno[,"SNP_ID"] == snp,"B_ALLELE"][1]
    geno[geno[,"SNP_ID"] == snp,"TRUE_GENO"] <- gsub("A", a.allele, geno[geno[,"SNP_ID"] == snp,"TRUE_GENO"])
    geno[geno[,"SNP_ID"] == snp,"TRUE_GENO"] <- gsub("B", b.allele, geno[geno[,"SNP_ID"] == snp,"TRUE_GENO"])
    geno[geno[,"SNP_ID"] == snp,"OBS_GENO"] <- gsub("A", a.allele, geno[geno[,"SNP_ID"] == snp,"OBS_GENO"])
    geno[geno[,"SNP_ID"] == snp,"OBS_GENO"] <- gsub("B", b.allele, geno[geno[,"SNP_ID"] == snp,"OBS_GENO"])
  }
  
  geno[,"INTENSITY_A"] <- cos(geno[,"OBS_ALLELIC_PROP"]*(pi/2))
  geno[,"INTENSITY_B"] <- sin(geno[,"OBS_ALLELIC_PROP"]*(pi/2))  
  
  return(geno[,c("SNP_ID", "B_ALLELE_FREQ", "TRUE_GENO", "SAMPLE_ID", "GENO_ERROR", "OBS_GENO", "OBS_ALLELIC_PROP", "INTENSITY_A", "INTENSITY_B", "A_ALLELE", "B_ALLELE")])
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

sim.geno.pool.fun <- function(sim.geno.indiv, map, indivs.in.sim.pools, n.in.pools) {
  
  #Get package
 # if("reshape2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("reshape2")} 
  library(reshape2)
  
  sel.indivs <- indivs.in.sim.pools[,"SAMPLE_ID"]
  
  #identify individuals in the pool
  #  if(length(unique(indivs.in.sim.pools[,"FAM_SET_ID"])) > 1) {
  #    sel.fam.agg <- sample(unique(indivs.in.sim.pools[,"FAM_SET_ID"]), n.in.pools, replace = FALSE, prob = NULL) 
  #  } else {
  #    sel.fam.agg <- unique(indivs.in.sim.pools[,"FAM_SET_ID"])
  #  }
  #  
  #  sel.indivs <- NULL
  #  for(i in sel.fam.agg) {
  #    x <- indivs.in.sim.pools[indivs.in.sim.pools[,"FAM_SET_ID"] == i,"SAMPLE_ID"]
  #    if(length(x) > 1) {
  #      sel.indiv <-  sample(x = x, size = 1)
  #    } else {
  #      sel.indiv <- x
  #    }
  #    rm(x)
  #    
  #    sel.indivs <- c(sel.indivs, sel.indiv)
  #   
  #  }
  #  rm(sel.indiv, sel.fam.agg)
  
  #generate geno.pool
  geno.pool <- sim.geno.indiv[sim.geno.indiv[,"SAMPLE_ID"] %in% sel.indivs,c("SNP_ID", "B_ALLELE_FREQ", "SAMPLE_ID", "TRUE_GENO")]
  geno.pool <- dcast(geno.pool, SNP_ID + B_ALLELE_FREQ ~ SAMPLE_ID, value.var = "TRUE_GENO") 
  
  #concatenate to get pooled genotypes
  if(n.in.pools > 1) {
    tmp <- c(geno.pool[,!colnames(geno.pool) %in% c("SNP_ID", "B_ALLELE_FREQ")], sep="")
    tmp <- do.call(paste, tmp)
    geno.pool <- as.data.frame(cbind(geno.pool[,c("SNP_ID", "B_ALLELE_FREQ")], tmp))
    rm(tmp)
  }
  
  colnames(geno.pool) <- c("SNP_ID", "B_ALLELE_FREQ","TRUE_ORDERED_GENO")
  
  #convert ordered genotypes to unordered
  ordered.genos <- expand.grid(rep(list(1:2), n.in.pools*2))
  ordered.genos[ordered.genos[,] == 1] <- "A"
  ordered.genos[ordered.genos[,] == 2] <- "B"
  tmp <- c(ordered.genos, sep="")
  ordered.genos <- do.call(paste, tmp)
  
  unordered.genos <- genotypes.fun(n = n.in.pools*2)
  for (i in 1:length(ordered.genos)) {
    unordered.genos[i] <- paste(sort(unlist(strsplit(ordered.genos[i], ""))), collapse = "")
  }
  unordered.genos <- as.data.frame(cbind(ordered.genos, unordered.genos))
  colnames(unordered.genos) <- c("TRUE_ORDERED_GENO","TRUE_UNORDERED_GENO")
  # geno.pool <- merge(geno.pool, unordered.genos, by = "TRUE_ORDERED_GENO", all.x = TRUE)
  geno.pool$TRUE_ORDERED_GENO    <- as.character(geno.pool$TRUE_ORDERED_GENO)
  unordered.genos$TRUE_ORDERED_GENO    <- as.character(unordered.genos$TRUE_ORDERED_GENO)
  geno.pool <- left_join(geno.pool, unordered.genos, by = "TRUE_ORDERED_GENO")
  geno.pool <- geno.pool[,c("SNP_ID", "B_ALLELE_FREQ","TRUE_UNORDERED_GENO")]
  colnames(geno.pool) <- c("SNP_ID", "B_ALLELE_FREQ", "TRUE_GENO")
  
  geno.pool[,"SAMPLE_ID"] <- max(sim.geno.indiv[,"SAMPLE_ID"]) + 1
  
  geno.pool <- geno.pool[order(geno.pool$SNP_ID),]
  
  #Use genotype replacement model of Marshall et al (1998)
  snp.id.error <-  map[runif(n = nrow(map)) < map[,"ERROR_RATE"],"SNP_ID"]  
  
  geno.pool[,"GENO_ERROR"] <- FALSE  
  geno.pool[,"OBS_GENO"]   <- geno.pool[,"TRUE_GENO"]
  
  if(length(snp.id.error) > 0) {
    
    replace <- NULL
    for(gene in 1:(n.in.pools*2)) {
      tmp <- runif(n = nrow(geno.pool)) < geno.pool[,"B_ALLELE_FREQ"]
      replace <- cbind(replace,tmp)
    }
    replace[replace[TRUE]] <- "B"
    replace[replace == "FALSE"] <- "A"
    replace <- as.data.frame(replace)
    colnames(replace) <- 1:ncol(replace)
    
    #concatenate to get pooled genotypes
    tmp <- c(replace, sep="")
    tmp <- do.call(paste, tmp)
    replace <- as.data.frame(cbind(geno.pool[,c("SNP_ID")], tmp))
    colnames(replace) <- c("SNP_ID", "OBS_ORDERED_GENO")
    rm(tmp)
    
    #convert ordered genotypes to unordered
    colnames(unordered.genos) <- c("OBS_ORDERED_GENO","OBS_UNORDERED_GENO")
    # replace <- merge(replace, unordered.genos, by = "OBS_ORDERED_GENO", all.x = TRUE)
    replace$OBS_ORDERED_GENO    <- as.character(replace$OBS_ORDERED_GENO)
    unordered.genos$OBS_ORDERED_GENO    <- as.character(unordered.genos$OBS_ORDERED_GENO)
    replace <- left_join(replace, unordered.genos, by = "OBS_ORDERED_GENO")
    
    replace <- replace[,c("SNP_ID","OBS_UNORDERED_GENO")]
    colnames(replace) <- c("SNP_ID", "OBS_GENO")
    replace <- replace[order(replace[,"SNP_ID"]),]
    
    snp.error <- geno.pool[,"SNP_ID"] %in% snp.id.error
    geno.pool[snp.error,"GENO_ERROR"] <- TRUE
    geno.pool <- geno.pool[order(geno.pool[,"SNP_ID"]),]
    geno.pool[snp.error, "OBS_GENO"]  <- replace[snp.error, "OBS_GENO"] 
    rm(replace, snp.error)
  }
  
  return(list(geno.pool = geno.pool, true.indivs.in.pool = sel.indivs))
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

sim.miss.in.snp.dat.indiv <- function(sim.snp.dat.indiv, map, missing.parents = NULL) {
  
  #map:
  #CHROMOSOME (integer), SNP_ID (character), GENETIC_POSITION (Morgans - numeric), 
  #  B_ALLELE_FREQ (numeric - 0-1), ERROR_RATE (numeric 0-1) , PROP_MISS (numeric 0-1)
  #  sim.snp.dat.indiv <- merge(sim.snp.dat.indiv, map[,c("SNP_ID", "PROP_MISS")], by = "SNP_ID", all.x = TRUE)
  sim.snp.dat.indiv$SNP_ID    <- as.character(sim.snp.dat.indiv$SNP_ID)
  
  map$SNP_ID    <- as.character(map$SNP_ID)
  sim.snp.dat.indiv <- left_join(sim.snp.dat.indiv, map[,c("SNP_ID", "PROP_MISS")], by = "SNP_ID")
  
  sim.snp.dat.indiv[,"MISSING"] <- runif(n = nrow(sim.snp.dat.indiv)) < sim.snp.dat.indiv[,"PROP_MISS"]
  
  #remove data from samples will all SNP missing
  if(!is.null(missing.parents)) {
    tmp <- as.character(sim.snp.dat.indiv[,"SAMPLE_ID"])
    missing.parents <- as.character(missing.parents)  
    sim.snp.dat.indiv[tmp %in% missing.parents,"MISSING"] <- TRUE
    rm(tmp)
  }
  
  sim.snp.dat.indiv[sim.snp.dat.indiv[,"MISSING"],colnames(sim.snp.dat.indiv) %in% c("INTENSITY_A",  "INTENSITY_B", "OBS_GENO", "GENOTYPE",  "OBS_ALLELIC_PROP")] <- NA
  #sim.snp.dat.indiv <- sim.snp.dat.indiv[!sim.snp.dat.indiv[,"MISSING"],] 
  sim.snp.dat.indiv <- sim.snp.dat.indiv[,!colnames(sim.snp.dat.indiv) %in% c("PROP_MISS", "MISSING")]  
  
  return(sim.snp.dat.indiv)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#' @export

sim.pool.fun <- function(map, ped, true.snp.param.indiv, indivs.in.sim.pools, missing.parents = NULL) {
  
  n.in.pools <- nrow(indivs.in.sim.pools)
  
  #Generate simulated data and manipulate data for max.likelihood.fun######################################
  
  #individual data
  geno.indiv <- sim.geno.indiv.fun(map = map, ped = ped)
  
  sim.snp.dat.indiv    <- sim.snp.dat.indiv.fun(geno = geno.indiv, true.snp.param.indiv = true.snp.param.indiv)
  sim.snp.dat.indiv    <- sim.miss.in.snp.dat.indiv(sim.snp.dat.indiv = sim.snp.dat.indiv, map, missing.parents = missing.parents)
  
  #pool data
  geno.pool    <- sim.geno.pool.fun(sim.geno.indiv = geno.indiv[,c("SAMPLE_ID", "B_ALLELE_FREQ", "SNP_ID", "TRUE_GENO")], 
                                    map = map,
                                    indivs.in.sim.pools = indivs.in.sim.pools, 
                                    n.in.pools = n.in.pools)
  true.indivs.in.pool  <- geno.pool$true.indivs.in.pool
  geno.pool     <- geno.pool$geno.pool  
  
  true.snp.param.pools <-  snp.param.pools.fun(snp.param.indiv = true.snp.param.indiv,
                                               n.in.pools) 
  sim.snp.dat.pools <- sim.snp.dat.indiv.fun(geno = geno.pool, true.snp.param.indiv = true.snp.param.pools)
  sim.snp.dat.pools  <- sim.miss.in.snp.dat.indiv(sim.snp.dat.indiv = sim.snp.dat.pools, map = map, missing.parents = missing.parents)  

  #estimate parameters from simulated data
  sim.snp.dat.indiv[,"GENOTYPE"] <- sim.snp.dat.indiv[,"OBS_GENO"]
  #  sim.snp.dat.indiv[sim.snp.dat.indiv[,"GENOTYPE"] == "AA" & !is.na(sim.snp.dat.indiv[,"GENOTYPE"]),"GENOTYPE"] <- "A"
  #  sim.snp.dat.indiv[sim.snp.dat.indiv[,"GENOTYPE"] == "CC" & !is.na(sim.snp.dat.indiv[,"GENOTYPE"]),"GENOTYPE"] <- "C"
  #  sim.snp.dat.indiv[sim.snp.dat.indiv[,"GENOTYPE"] == "GG" & !is.na(sim.snp.dat.indiv[,"GENOTYPE"]),"GENOTYPE"] <- "G"
  #  sim.snp.dat.indiv[sim.snp.dat.indiv[,"GENOTYPE"] == "TT" & !is.na(sim.snp.dat.indiv[,"GENOTYPE"]),"GENOTYPE"] <- "T"
  
  return(list(geno.indiv = geno.indiv,
              geno.pool = geno.pool,
              true.indivs.in.pool = true.indivs.in.pool,
              sim.snp.dat.indiv = sim.snp.dat.indiv,
              sim.snp.dat.pools = sim.snp.dat.pools,
              true.snp.param.pools = true.snp.param.pools
  ))
  
}






