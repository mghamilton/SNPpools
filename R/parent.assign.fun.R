# June 2020
# Matthew Hamilton

#' parent.assign.fun
#' 
#' This function assigns parents to pooled samples using one of seven approaches:
#' \itemize{
#'  \item{'Least squares' method outlined in Henshall et al. 2014 (method = "Least_squares")} 
#'  \item{'Least squares minimum sum of squares' method whereby the sum of squares of all parental combinations are computed and 
#'   the combination with the minimum value identified.  Refer to Hamilton 2020 (method = "Least_squares", beta.min.ss = TRUE)} 
#'  \item{'Quantitative maximum likelihood' method whereby the approach using quantititive genotypes for parentatge assignment
#'  outlined in Henshall et al. 2014 is extended to pooled DNA samples.  Refer to Hamilton 2020 (method = "Quantitative")}
#'  \item{'Discrete maximum likelihood from genotype probabilities' method whereby discrete genotypes are derived from genotype probabilities and  
#'  used for parentatge assignment of pooled DNA samples by extending the discrete genotype maximum likelihood approach outlined in 
#'  Henshall et al. 2014 to pooled DNA samples.  Refer to Hamilton 2020 (method = "Discrete", discrete.method = "geno.probs").}
#'  \item{'Discrete maximum likelihood from genotype assignments' method whereby discrete genotypes are provided as input and  
#'  used for parentatge assignment of pooled DNA samples by extending the discrete genotype maximum likelihood approach outlined in 
#'  Henshall et al. 2014 to pooled DNA samples.  Refer to Hamilton 2020 (method = "Discrete", discrete.method = "assigned.genos")}
#'  \item{'Exclusion from genotype probabilities' method whereby discrete genotypes are derived from genotype probabilities and  
#'  used for parentatge assignment of pooled DNA samples by extending the exclusion approach outlined in 
#'  Henshall et al. 2014 to pooled DNA samples  Refer to Hamilton 2020 (method = "Exclusion", discrete.method = "geno.probs").}
#'  \item{'Exclusion from genotype assignments' method whereby discrete genotypes are provided as input and  
#'  used for parentatge assignment of pooled DNA samples by extending the exclusion approach outlined in 
#'  Henshall et al. 2014 to pooled DNA samples.  Refer to Hamilton 2020 (method = "Exclusion", discrete.method = "assigned.genos")}
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
#' @param snp.dat.pools is a data frame with the following headings.  Note that all pooled DNA samples 
#' in this dataframe must be comprised of DNA from the same number of individuals (see n.in.pools) (class in parentheses):
#'  \itemize{
#'  \item{SAMPLE_ID is the pooled sample identifier (integer).}
#'  \item{SNP_ID is the SNP identifier (character).}
#'  \item{INTENSITY_A is the signal intensity for allele A. Not required if method does not include 'Quantitative' or 'Least_squares' and discrete.method  = "geno.probs" (numeric).}
#'  \item{INTENSITY_B is the signal intensity for allele B. Not required if method does not include 'Quantitative' or 'Least_squares' and discrete.method  = "geno.probs" (numeric).}
#'  \item{GENOTYPE is the assigned unordered genotype. Not required if discrete.method = "geno.probs" (character).}
#' } 
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
#' @param snp.error.underlying is not used if snp.error.assumed is not NULL (default = NULL). Must be either:
#' \itemize{
#'  \item{NULL.}
#'  \item{a numeric variable between 0 and 1 inclusive.  Used to comptute SNP_ERROR_TILDE from SNP_ERROR_HAT according
#'                      to the approach outlined on the left of page 5 of Henshall et al. 2014 using individual 
#'                      (i.e. not pooled) data only.  If snp.error.underlying = 0 then SNP_ERROR_TILDE = SNP_ERROR_HAT.}
#' } 
#' @param snp.param.indiv is the output of snp.param.indiv.fun (default = NULL).  That is, it is a data frame with the following 
#' headings (class in parentheses):
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
#' @param snp.param.pools is the output of snp.param.pools.fun.  That is, it is a data frame with the following headings (class in parentheses)
#' \itemize{
#'  \item{'SNP_ID' is the SNP identifier (character).} 
#'  \item{'MEAN_P_AAAA' is the mean of allelic proportion for homozygous A genotypes (numeric).}
#'  \item{'SD_P_AAAA' is the standard deviation of allelic proportion for homozygous A genotypes (numeric).} 
#'  \item{'MEAN_P_AAAB' is the mean of allelic proportion for unordered AAAB genotypes (numeric).}
#'  \item{'SD_P_AAAB' is the standard deviation of allelic proportion for unordered AAAB genotypes (numeric).} 
#'  \item{'MEAN_P_AABB' is the mean of allelic proportion for unordered AABB genotypes (numeric).}
#'  \item{'SD_P_AABB' is the standard deviation of allelic proportion for unordered AABB genotypes (numeric).}  
#'  \item{'MEAN_P_ABBB' is the mean of allelic proportion for unordered ABBB genotypes (numeric).}
#'  \item{'SD_P_ABBB' is the standard deviation of allelic proportion for unordered ABBB genotypes (numeric).}
#'  \item{'MEAN_P_BBBB' is the mean of allelic proportion for homozygous B genotypes (numeric).}
#'  \item{'SD_P_BBBB' is the standard deviation of allelic proportion for homozygous B genotypes (numeric).}   
#'  \item{'A_ALLELE' is the base represented by allele A (i.e. 'A', 'C', 'G' or 'T') (character).}
#'  \item{'B_ALLELE' is the base represented by allele B (i.e. 'A', 'C', 'G' or 'T') (character).}
#' } 
#' @param min.sd is a numberic variable defining a lower bound to be applied to estimates of the 
#' standard deviation of allelic proportion for genotypes in snp.param.indiv and snp.param.pools (default = 0)
#' @param fams is a data frame with the following headings (class in parentheses):
#' \itemize{
#'  \item{'FAMILY_ID' is the family identifier (integer).} 
#'  \item{'SIRE_ID' is the sire identifier (integer).} 
#'  \item{'DAM_ID' is the dam identifier (integer).} 
#' } 
#' @param fam.set.combns is a data frame with the following headings (class in parentheses). Note: if fam.set.combns = NULL
#' (see 'pooling by phenotype' example below), 
#' FAMILY_ID is taken from the 'fams' and duplicated n.in.pools times, FAM_SET_ID = 1 for the first duplication of 
#' FAMILY_IDs, 2 for the second etc and FAM_SET_COMBN_ID = 1 (default = NULL):
#' \itemize{
#'  \item{'FAM_SET_COMBN_ID' is the family set combination identifier (integer). A 'family set combination' is a combination of 'family sets'.  
#'  Each pooled sample must be associated with one only family set combination but a family set combination
#'  may be assoicated with multiple pooled samples using the fam.set.combns.by.pool input below.} 
#'  \item{'FAM_SET_ID' is the family set identifier (integer).  A 'family set' is a group of families of which one is known to be the true family 
#'  of one of the individuals in a pooled sample.  Within each 'family set combination' there must be a 'family set'
#'  for each individual in a pooled sample (i.e. if n.in.pools = 2 there must be two family sets in each family set combination)} 
#'  \item{'FAMILY_ID' is the family identifier (integer).} 
#' } 
#' @param fam.set.combns.by.pool is a data frame linking pooled samples with family set combinations.  
#' It has the following headings (class in parentheses).  Note: if fam.set.combns is NULL (see 'pooling by phenotype' example below),
#' fam.set.combns.by.pool is made NULL.
#' If fam.set.combns.by.pool = NULL, FAM_SET_COMBN_ID = 1 and SAMPLE_ID is taken from the 
#' 'snp.dat.pools' input (default = NULL):
#' \itemize{
#'  \item{SAMPLE_ID is the pooled sample identifier (integer).}
#'  \item{'FAM_SET_COMBN_ID' is the family set combination identifier (integer).} 
#' } 
#' @param skip.checks is a logical variable.  If FALSE parent.assign.fun data checks are not undertaken. 
#' @return  
#' \cr
#' \strong{Primary outputs}
#' @return \code{most.like.parents.quant}  
#' \cr \cr
#' Applicable when method = "Quantitative".  Identifies the most likely parental combination and delta LODs for 
#' individual parents.  Second most likely (alternative) parents are also presented.  Refer to Hamilton 2020.  Example fields for n.in.pools = 2:
#' \itemize{
#'  \item{SAMPLE_ID (integer).}
#'  \item{PARENT_COMBN_ID (integer).}
#'  \item{MISS_PARENT_SNP_DATA_PROP (numeric).}
#'  \item{MISS_POOL_SNP_DATA_PROP (numeric).}
#'  \item{NO_MISS_PARENT_OR_POOL_PROP (numeric).}
#'  \item{MIN_LOGL (numeric).}
#'  \item{MIN_LOGL_SNP (character).}
#'  \item{MAX_LOGL (numeric).}
#'  \item{MAX_LOGL_SNP (character).}
#'  \item{RANGE_5_TO_95_LOGL (numeric).}
#'  \item{LOD (numeric).}
#'  \item{FAM_SET_COMBN_ID (integer).}
#'  \item{PARENT_1 (integer).}
#'  \item{PARENT_2 (integer).}
#'  \item{PARENT_3 (integer).}
#'  \item{PARENT_4 (integer).}
#'  \item{FAM_COMBN_ID (integer).}
#'  \item{FAMILY_ID_1 (integer).}
#'  \item{FAMILY_ID_2 (integer).}
#'  \item{PARENT_1_DELTA_LOD (logical).}
#'  \item{PARENT_2_DELTA_LOD (logical).}
#'  \item{PARENT_3_DELTA_LOD (numeric).}
#'  \item{PARENT_4_DELTA_LOD (numeric).}
#'  \item{ALT_PARENT_1 (logical).}
#'  \item{ALT_PARENT_2 (logical).}
#'  \item{ALT_PARENT_3 (numeric).}
#'  \item{ALT_PARENT_4 (numeric).}
#'  \item{ALT_PARENT_COMBN_1 (logical).}
#'  \item{ALT_PARENT_COMBN_2 (logical).}
#'  \item{ALT_PARENT_COMBN_3 (integer).}
#'  \item{ALT_PARENT_COMBN_4 (integer).}
#'  \item{ALT_FAM_COMBN_1 (integer).}
#'  \item{ALT_FAM_COMBN_2 (integer).}
#'  \item{ALT_FAM_COMBN_3 (integer).}
#'  \item{ALT_FAM_COMBN_4 (integer).}
#' }
#' \cr
#' @return \code{most.like.parents.discrete}
#' \cr  \cr
#' Applicable when method = "Discrete".  Identifies the most likely parental 
#' combination and delta LODs for individual parents.  Second most likely (alternative) parents are also presented.  
#' Refer to Hamilton 2020.  Example fields for n.in.pools = 2:
#' \itemize{
#'  \item{SAMPLE_ID (integer).}
#'  \item{PARENT_COMBN_ID (integer).}
#'  \item{MISS_PARENT_SNP_DATA_PROP (numeric).}
#'  \item{MISS_POOL_SNP_DATA_PROP (numeric).}
#'  \item{NO_MISS_PARENT_OR_POOL_PROP (numeric).}
#'  \item{MIN_LOGL (numeric).}
#'  \item{MIN_LOGL_SNP (character).}
#'  \item{MAX_LOGL (numeric).}
#'  \item{MAX_LOGL_SNP (character).}
#'  \item{RANGE_5_TO_95_LOGL (numeric).}
#'  \item{LOD (numeric).}
#'  \item{FAM_SET_COMBN_ID (integer).}
#'  \item{PARENT_1 (integer).}
#'  \item{PARENT_2 (integer).}
#'  \item{PARENT_3 (integer).}
#'  \item{PARENT_4 (integer).}
#'  \item{FAM_COMBN_ID (integer).}
#'  \item{FAMILY_ID_1 (integer).}
#'  \item{FAMILY_ID_2 (integer).}
#'  \item{PARENT_1_DELTA_LOD (logical).}
#'  \item{PARENT_2_DELTA_LOD (logical).}
#'  \item{PARENT_3_DELTA_LOD (numeric).}
#'  \item{PARENT_4_DELTA_LOD (numeric).}
#'  \item{ALT_PARENT_1 (logical).}
#'  \item{ALT_PARENT_2 (logical).}
#'  \item{ALT_PARENT_3 (numeric).}
#'  \item{ALT_PARENT_4 (numeric).}
#'  \item{ALT_PARENT_COMBN_1 (logical).}
#'  \item{ALT_PARENT_COMBN_2 (logical).}
#'  \item{ALT_PARENT_COMBN_3 (integer).}
#'  \item{ALT_PARENT_COMBN_4 (integer).}
#'  \item{ALT_FAM_COMBN_1 (integer).}
#'  \item{ALT_FAM_COMBN_2 (integer).}
#'  \item{ALT_FAM_COMBN_3 (integer).}
#'  \item{ALT_FAM_COMBN_4 (integer).}
#' }
#' \cr
#' @return \code{most.like.parents.excl}  
#' \cr  \cr
#' Applicable when method = "Exclusion".  Identifies the most likely parental 
#' combination.  Refer to Hamilton 2020.  Second most likely (alternative) parental 
#' combination is also presented.  
#' Example fields for n.in.pools = 2:
#' \itemize{
#'  \item{SAMPLE_ID (integer).}
#'  \item{PARENT_COMBN_ID (integer).}
#'  \item{MISMATCHES (integer).}
#'  \item{SNP_COUNT (integer).}
#'  \item{MISMATCH_PROP (numeric).}
#'  \item{MISMATCH_PROP_SE (numeric).}
#'  \item{MISMATCH_PROP_Z (numeric).}
#'  \item{FAM_COMBN_ID (integer).}
#'  \item{FAMILY_ID_1 (integer).}
#'  \item{FAMILY_ID_2 (integer).}
#'  \item{PARENT_1 (integer).}
#'  \item{PARENT_2 (integer).}
#'  \item{PARENT_3 (integer).}
#'  \item{PARENT_4 (integer).}
#'  \item{ALT_PARENT_COMBN_ID (integer).}
#'  \item{ALT_FAM_COMBN_ID (integer).}
#'  \item{ALT_MISMATCHES (integer).}
#'  \item{ALT_SNP_COUNT (integer).}
#'  \item{ALT_MISMATCH_PROP (numeric).}
#'  \item{ALT_MISMATCH_PROP_SE (numeric).}
#'  \item{ALT_MISMATCH_PROP_Z (numeric).}
#' }
#' \cr
#' @return \code{most.like.parents.excl.non.dup}
#' \cr \cr
#'  Applicable when method = "Exclusion".  Identifies the most likely parental
#'  combination - simplified output with multiple combinations with the same number of mismatches (duplicated 
#'  SAMPLE_IDs) removed.  Refer to Hamilton 2020.  Example fields for n.in.pools = 2:
#' \itemize{
#'  \item{SAMPLE_ID (integer).}
#'  \item{MISMATCHES (integer).}
#'  \item{SNP_COUNT (integer).}
#'  \item{MISMATCH_PROP (numeric).}
#'  \item{MISMATCH_PROP_SE (numeric).}
#'  \item{MISMATCH_PROP_Z (numeric).}
#'  \item{FAMILY_ID_1 (integer).}
#'  \item{FAMILY_ID_2 (integer).}
#'  \item{PARENT_1 (integer).}
#'  \item{PARENT_2 (integer).}
#'  \item{PARENT_3 (integer).}
#'  \item{PARENT_4 (integer).}
#' }
#' \cr
#' @return \code{beta}
#' \cr \cr
#' Applicable when method = "Least_squares". Identifies the most likely parental combination:
#' \itemize{
#'  \item{SAMPLE_ID (integer).}
#'  \item{SIRE_ID (integer).}
#'  \item{DAM_ID (integer).}
#'  \item{FAMILY_ID (integer).}
#'  \item{BETA_STAR Refer to Henshall et al. 2014  (numeric).}
#'  \item{BETA_HAT Refer to Henshall et al. 2014 (numeric).}
#'  \item{BETA_HAT_CONSTRAINED  Constrained to have equal contributions from each FAMILY_SET_ID. Refer to Hamilton 2020}
#'  \item{BETA_MIN_SS Applicable when beta.min.ss = "TRUE". Constrained beta with minimum sum of squares is retained.  Beta constrained to have equal contributions from each FAMILY_SET_ID.  Refer to Hamilton 2020}
#' }
#' \cr
#' \strong{Primary plots}
#' @return \code{bar.png} 
#' \itemize{
#'  \item{Bar plot of BETA_HAT by FAMILY_ID.  Applicable to 'Least_squares' method only.  Output in a directory named 'Results' on the current working directory.} 
#' }
#' @return \code{discrete.png} 
#' \itemize{
#'  \item{Scatter plot of '5-95 percentile range of log-likelihood ratios' against 'Log odds (LOD) scores' for each 
#'  possible family combination.  Ideally there is an isolated point in the bottom right represtenting the 
#'  correct family combination.  Applicable to 'Discrete' method only.  Output in a directory named 
#'  'Results/lod.scatter' on the current working directory.} 
#' }
#' @return \code{quantitative.png} 
#' \itemize{
#'  \item{Scatter plot of '5-95 percentile range of log-likelihood ratios' against 'Log odds (LOD) scores' for each 
#'  possible family combination.  Ideally there is an isolated point in the bottom right represtenting the 
#'  correct family combination.  Applicable to 'Quantitative' method only.  Output in a directory named 
#'  'Results/lod.scatter' on the current working directory.} 
#' }
#' \cr
#' \strong{Intermediate outputs}
#' \cr
#' @return 'Dij' Applicable when method = "Discrete" or "Exclusion". Refer to Henshall et al. 2014
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AA_GENO_PROB (numeric).}
#  \item{AB_GENO_PROB (numeric).}
#  \item{BA_GENO_PROB (numeric).}
#  \item{BB_GENO_PROB (numeric).}
#  \item{A_TRANS_PROB (numeric).}
#  \item{B_TRANS_PROB (numeric).}
# }
#' @return 'dkj' Applicable when method = "Discrete" or "Exclusion". Refer to Hamilton 2020
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
# }
#' @return 'dklj.adj' Applicable when method = "Discrete" or "Exclusion". Refer to Hamilton 2020 (dkj.star)
#' #:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{SAMPLE_ID (integer).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
# }
#' @return 'fkj.and.weight' Applicable when method = "Least_squares". Refer to Henshall et al. 2014
#' #:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{SAMPLE_ID (integer).}
#  \item{MEAN_P_AA (numeric).}
#  \item{MEAN_P_AB (numeric).}
#  \item{MEAN_P_BB (numeric).}
#  \item{ALLELIC_PROP_POOL (numeric).}
#  \item{FREQ_POOL (numeric).}
#  \item{FREQ_POOL_ERROR_WT (numeric).}
# }
#' @return 'Gij' Applicable when method = "Discrete" or ""Quantitative". Refer to Henshall et al. 2014
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AA_GENO_PROB (numeric).}
#  \item{AB_GENO_PROB (numeric).}
#  \item{BA_GENO_PROB (numeric).}
#  \item{BB_GENO_PROB (numeric).}
#  \item{A_TRANS_PROB (numeric).}
#  \item{B_TRANS_PROB (numeric).}
#  \item{ALLELIC_PROP_INDIV (numeric).}
# }
#' @return 'gkj' Applicable when method = "Discrete" or ""Quantitative". Refer to Hamilton 2020 
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
# }
#' @return 'gklj.adj' Applicable when method = "Discrete" or ""Quantitative". Refer to Hamilton 2020 (gkj.star)
#' #:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{SAMPLE_ID (integer).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
# }
#' @return 'flj.probs' Applicable when discrete.method = "geno.probs". Refer to Hamilton 2020 (fj)
#' #:
# \itemize{
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{A_TRANS_PROB (numeric).}
#  \item{B_TRANS_PROB (numeric).}
#  \item{SAMPLE_ID (integer).}
# }
#' @return 'flj.geno' Applicable when discrete.method = "assigned.genos". Refer to Hamilton 2020 (fj)
#' #:
# \itemize{
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{A_TRANS_PROB (numeric).}
#  \item{B_TRANS_PROB (numeric).}
#  \item{SAMPLE_ID (integer).}
# }
#' @return 'lambda.kj' Applicable when method = "Discrete" or ""Quantitative". Refer to Hamilton 2020. 
#' #Example output for n.in.pools = 2#:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AAAA_LAMBDA (numeric).}
#  \item{AAAB_LAMBDA (numeric).}
#  \item{AABB_LAMBDA (numeric).}
#  \item{ABBB_LAMBDA (numeric).}
#  \item{BBBB_LAMBDA (numeric).}
#  \item{ALLELIC_PROP_POOL (numeric).}
# }
#' @return 'lod.duos.discrete' Applicable when method = "Discrete". Refer to Hamilton 2020.
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{MISS_PARENT_SNP_DATA_PROP (numeric).}
#  \item{MISS_POOL_SNP_DATA_PROP (numeric).}
#  \item{NO_MISS_PARENT_OR_POOL_PROP (numeric).}
#  \item{MIN_LOGL (numeric).}
#  \item{MIN_LOGL_SNP (character).}
#  \item{MAX_LOGL (numeric).}
#  \item{MAX_LOGL_SNP (character).}
#  \item{RANGE_5_TO_95_LOGL (numeric).}
#  \item{LOD (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{PARENT_1 (integer).}
#  \item{PARENT_2 (integer).}
#  \item{PARENT_3 (integer).}
#  \item{PARENT_4 (integer).}
# }
#' @return 'lod.duos.quant' Applicable when method = "Quantitative". Refer to Hamilton 2020.
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{MISS_PARENT_SNP_DATA_PROP (numeric).}
#  \item{MISS_POOL_SNP_DATA_PROP (numeric).}
#  \item{NO_MISS_PARENT_OR_POOL_PROP (numeric).}
#  \item{MIN_LOGL (numeric).}
#  \item{MIN_LOGL_SNP (character).}
#  \item{MAX_LOGL (numeric).}
#  \item{MAX_LOGL_SNP (character).}
#  \item{RANGE_5_TO_95_LOGL (numeric).}
#  \item{LOD (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{PARENT_1 (integer).}
#  \item{PARENT_2 (integer).}
#  \item{PARENT_3 (integer).}
#  \item{PARENT_4 (integer).}
# }
#' @return 'logl.duos.discrete' Applicable when method = "Discrete". Refer to Hamilton 2020. Only outputted for the final SAMPLE_ID
#' #:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{SAMPLE_ID (integer).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{LIKE_(numeric)ERATOR (numeric).}
#  \item{LIKE_DENOMINATOR (numeric).}
#  \item{LOG_LIKE_RATIO (numeric).}
#  \item{MISS_PARENT_SNP_DATA_PROP (numeric).}
#  \item{MISS_POOL_SNP_DATA_PROP (numeric).}
# }
#' @return 'logl.duos.quant' Applicable when method = "Quantitative". Refer to Hamilton 2020. Only outputted for the final SAMPLE_ID
#' #:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{SAMPLE_ID (integer).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{LIKE_(numeric)ERATOR (numeric).}
#  \item{LIKE_DENOMINATOR (numeric).}
#  \item{LOG_LIKE_RATIO (numeric).}
#  \item{MISS_PARENT_SNP_DATA_PROP (numeric).}
#  \item{MISS_POOL_SNP_DATA_PROP (numeric).}
# }
#' @return 'mismatches' Applicable when method = "Exclusion". Refer to Hamilton 2020.
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{MISMATCHES (integer).}
#  \item{SNP_COUNT (integer).}
#  \item{MISMATCH_PROP (numeric).}
#  \item{MISMATCH_PROP_SE (numeric).}
#  \item{MISMATCH_PROP_Z (numeric).}
# }
#' @return 'mismatches.by.snp' Applicable when method = "Exclusion". Refer to Hamilton 2020. 
#' #Example output for n.in.pools = 2:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{AAAA_TRANS (numeric).}
#  \item{AAAB_TRANS (numeric).}
#  \item{AABB_TRANS (numeric).}
#  \item{ABBB_TRANS (numeric).}
#  \item{BBBB_TRANS (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{MISS_PARENT_COUNT (numeric).}
#  \item{MISS_PARENT (logical).}
#  \item{SAMPLE_ID (integer).}
#  \item{AAAA_POOLS (numeric).}
#  \item{AAAB_POOLS (numeric).}
#  \item{AABB_POOLS (numeric).}
#  \item{ABBB_POOLS (numeric).}
#  \item{BBBB_POOLS (numeric).}
#  \item{MISS_POOL (logical).}
#  \item{MISMATCHES (integer).}
# }
#' @return 'nlj.probs' Applicable when discrete.method = "geno.probs". Refer to Hamilton 2020 (nj). 
#' #Example output for n.in.pools = 2:
# \itemize{
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
# }
#' @return 'nlj.geno'  Applicable when discrete.method = "geno.probs". Refer to Hamilton 2020 (nj).
#' # Example output for n.in.pools = 2:
# \itemize{
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
# }
#' @return 'parent.combns':
# \itemize{
#  \item{FAM_COMBN_ID (integer).}
#  \item{FAMILY_ID_1 (integer).}
#  \item{FAMILY_ID_2 (integer).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{PARENT_1 (integer).}
#  \item{PARENT_2 (integer).}
#  \item{PARENT_3 (integer).}
#  \item{PARENT_4 (integer).}
# }
#' @return 'phi.ij' Applicable when method = "Discrete" or ""Quantitative". Refer to Henshall et al. 2014
#' #:
# \itemize{
#  \item{SAMPLE_ID (integer).}
#  \item{SNP_ID (character).}
#  \item{AA_PHI (numeric).}
#  \item{AB_PHI (numeric).}
#  \item{BA_PHI (numeric).}
#  \item{BB_PHI (numeric).}
#  \item{ALLELIC_PROP_INDIV (numeric).}
# }
#' @return 'snp.error.probs' Applicable when discrete.method = "geno.probs". Refer to Henshall et al. 2014
#' #:
# \itemize{
#  \item{SNP_ID Factor}
#  \item{SNP_ERROR_TILDE (numeric).}
# }
#' @return 'snp.error.geno'  Applicable when discrete.method = "assigned.genos". Refer to Henshall et al. 2014
#' #:
# \itemize{
#  \item{SNP_ID Factor}
#  \item{SNP_ERROR_TILDE (numeric).}
# }
#' @return 'tclj.adj.quant' Applicable when method = "Quantitative". Refer to Hamilton 2020 (tcj.star). 
#' #Example output for n.in.pools = 2:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
# }
#' @return 'tclj.adj.discrete' Applicable when method = "Discrete". Refer to Hamilton 2020 (tcj.star). 
#' #Example output for n.in.pools = 2:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
# }
#' @return 'tclj.discrete' Applicable when method = "Discrete". Refer to Hamilton 2020 (tcj). 
#' #Example output for n.in.pools = 2:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
# }
#' @return 'tclj.ls' Applicable when method = "Least_squares". Refer to Hamilton 2020 (tcj).
#' 3:
# \itemize{
#  \item{SNP_ID (character).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{MISS_PARENT_COUNT (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
#  \item{B_TRANS_PROB (numeric).}
# }
#' @return 'tclj.quant' Applicable when method = "Quantitative". Refer to Hamilton 2020 (tcj). 
#' #Example output for n.in.pools = 2::
# \itemize{
#  \item{SNP_ID (character).}
#  \item{PARENT_COMBN_ID (integer).}
#  \item{AAAA (numeric).}
#  \item{AAAB (numeric).}
#  \item{AABB (numeric).}
#  \item{ABBB (numeric).}
#  \item{BBBB (numeric).}
#  \item{FAM_SET_COMBN_ID (integer).}
# }
#' @return 'Xl.mat' List Applicable when method = "Quantitative". Refer to Henshall et al 2014 (X)
#' #:
# \itemize{
#  \item{ (numeric).}
# }

#' @examples
#' 
#' #' #Retrieve data for 'pooling by phenotype' example from Hamilton 2020
#' data(shrimp.snp.dat.indiv)
#' data(shrimp.snp.dat.pools)
#' data(shrimp.fams)
#' 
#' #Compute SNP parameters
#' shrimp.snp.param.indiv <- snp.param.indiv.fun(shrimp.snp.dat.indiv)
#' shrimp.snp.param.pools <- snp.param.pools.fun(shrimp.snp.param.indiv, n.in.pools = 2)
#'
#' #Assign parentage using the quantitative maximum likelihood method
#' parent.assign.fun(method= "Quantitative",
#'                   snp.dat.indiv = shrimp.snp.dat.indiv, 
#'                   snp.dat.pools = shrimp.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = shrimp.snp.param.indiv,
#'                   snp.param.pools = shrimp.snp.param.pools,                  
#'                   fams = shrimp.fams)  
#' 
#' #Retrieve data for 'pooling for individual parentage assignment' example from Hamilton 2020
#' data(ab.snp.dat.indiv)
#' data(ab.snp.dat.pools)
#' data(ab.fams)
#' data(ab.fam.set.combns)
#' data(ab.fam.set.combns.by.pool)
#' 
#' #Compute SNP parameters
#' ab.snp.param.indiv <- snp.param.indiv.fun(ab.snp.dat.indiv)
#' ab.snp.param.pools <- snp.param.pools.fun(ab.snp.param.indiv, n.in.pools = 3)
#'
#' #Assign parentage using the quantitative maximum likelihood method
#' parent.assign.fun(method= "Quantitative",
#'                   snp.dat.indiv = ab.snp.dat.indiv, 
#'                   snp.dat.pools = ab.snp.dat.pools,
#'                   n.in.pools = 3,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = ab.snp.param.indiv,
#'                   snp.param.pools = ab.snp.param.pools,                  
#'                   fams = ab.fams,
#'                   fam.set.combns = ab.fam.set.combns,
#'                   fam.set.combns.by.pool = ab.fam.set.combns.by.pool) 
#'                   
#' #Retrieve data for small worked example from Hamilton 2020
#' data(Ham.snp.dat.indiv)
#' data(Ham.snp.dat.pools)
#' data(Ham.fams)
#' data(Ham.fam.set.combns)
#' data(Ham.fam.set.combns.by.pool)
#' 
#' #Compute SNP parameters
#' Ham.snp.param.indiv <- snp.param.indiv.fun(Ham.snp.dat.indiv)
#' Ham.snp.param.pools <- snp.param.pools.fun(Ham.snp.param.indiv, n.in.pools = 2)
#' 
#' #Assign parentage using the least squares method
#' parent.assign.fun(method = "Least_squares",
#'                   beta.min.ss = TRUE, 
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = Ham.snp.param.indiv,
#'                   snp.param.pools = Ham.snp.param.pools,                  
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)
#'                   
#' #Assign parentage using the quantitative maximum likelihood method
#' parent.assign.fun(method= "Quantitative",
#'                   threshold.indiv = 0.98,         
#'                   threshold.pools = 0.98,         
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = Ham.snp.param.indiv,
#'                   snp.param.pools = Ham.snp.param.pools,                  
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)  
#'                   
#' #Assign parentage using the discrete maximum likelihood method 
#' #(discrete.method = "geno.probs")
#' parent.assign.fun(method= "Discrete",
#'                   discrete.method = "geno.probs",
#'                   threshold.indiv = 0.98,         
#'                   threshold.pools = 0.98,         
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = Ham.snp.param.indiv,
#'                   snp.param.pools = Ham.snp.param.pools,                  
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)  
#'                   
#' #Assign parentage using the discrete maximum likelihood method 
#' #(discrete.method = "assigned.genos")
#' parent.assign.fun(method= "Discrete",
#'                   discrete.method = "assigned.genos",
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)    
#'   
#' #Assign parentage using the exclusion method 
#' #(discrete.method = "geno.probs")
#' parent.assign.fun(method= "Exclusion",
#'                   discrete.method = "geno.probs",
#'                   threshold.indiv = 0.98,         
#'                   threshold.pools = 0.98,         
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = Ham.snp.param.indiv,
#'                   snp.param.pools = Ham.snp.param.pools,                  
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)   
#'                                   
#' #Assign parentage using the exclusion method 
#' #(discrete.method = "assigned.genos")
#' parent.assign.fun(method= "Exclusion",
#'                   discrete.method = "assigned.genos",
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)   
#'  
#' #Assign parentage using multiple methods
#' #(discrete.method = "geno.probs")
#' parent.assign.fun(method = c("Least_squares", "Quantitative", "Discrete", "Exclusion"),
#'                   beta.min.ss = TRUE, 
#'                   discrete.method = "geno.probs",
#'                   threshold.indiv = 0.98,         
#'                   threshold.pools = 0.98,         
#'                   snp.dat.indiv = Ham.snp.dat.indiv, 
#'                   snp.dat.pools = Ham.snp.dat.pools,
#'                   n.in.pools = 2,
#'                   snp.error.assumed = 0.01,
#'                   snp.param.indiv = Ham.snp.param.indiv,
#'                   snp.param.pools = Ham.snp.param.pools,                  
#'                   fams = Ham.fams,
#'                   fam.set.combns = Ham.fam.set.combns,
#'                   fam.set.combns.by.pool = Ham.fam.set.combns.by.pool)  
#'                   
#' @references Henshall JM, Dierens, L Sellars MJ (2014) Quantitative analysis of low-density SNP data for parentage assignment and estimation of family contributions to pooled samples. Genetics Selection Evolution 46, 51. https://doi 10.1186/s12711-014-0051-y 
#' @references Hamilton MG (2020) Maximum likelihood parentage assignment using quantitative genotypes

#' @export

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

parent.assign.fun <- function(method, #c("Quantitative", "Discrete", "Exclusion", "Least_squares"),
                              beta.min.ss = FALSE, #Appicable to least_squares method.  If TRUE, beta constrained to integers according to n.in.pools and minimum sum of squares identified
                              discrete.method = "geno.probs", #"geno.probs" or "assigned.genos" 
                              threshold.indiv = NULL,         #Dij.from.Gij.fun.  Not used if discrete.method = "assigned.genos"
                              threshold.pools = NULL,         #dkj.from.gkj.fun.  Not used if discrete.method = "assigned.genos"
                              
                              #SNP data
                              snp.dat.indiv, 
                              snp.dat.pools,
                              n.in.pools,
                              min.intensity        = 0,    #pij.fun.  
                              snp.error.assumed    = NULL, #If not null then this error is applied to all SNP.
                              snp.error.underlying = NULL, #adj.geno.prob.fun.  Not required if snp.error.assumed is not NULL.
                              
                              #SNP parameters
                              snp.param.indiv = NULL,
                              snp.param.pools = NULL,
                              min.sd          = 0,         #phi.ij.fun   lambda.ij.fun  
                              
                              #Define pools
                              fams,
                              fam.set.combns = NULL,
                              fam.set.combns.by.pool = NULL,

                             skip.checks = FALSE
                             ) {

  print(Sys.time())
  print("Running parent.assign.fun")  
  
  # load required packages
  if("dplyr" %in% installed.packages()[, "Package"] == FALSE) {install.packages("dplyr", repos='https://cran.csiro.au/')} 
  library(dplyr) 
  if("mgcv" %in% installed.packages()[, "Package"] == FALSE) {install.packages("mgcv", repos='https://cran.csiro.au/')} 
  library(mgcv)
  if("RColorBrewer" %in% installed.packages()[, "Package"] == FALSE) {install.packages("RColorBrewer", repos='https://cran.csiro.au/')} 
  library(RColorBrewer)
  if("gplots" %in% installed.packages()[, "Package"] == FALSE) {install.packages("gplots", repos='https://cran.csiro.au/')} 
  library(gplots)
  if("ggplot2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("ggplot2", repos='https://cran.csiro.au/')} 
  library(ggplot2)
  if("reshape2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("reshape2", repos='https://cran.csiro.au/')} 
  library(reshape2)
  
  #Define plot parameters for bar.plot.fun - these were inputs in the parent.assign.fun but have been moved to simplify inputs
  file.name = ""              #Text. Name of thebar plot file.  Only relevant if method includes 'Least_squares'
  var = "BETA_HAT"            #Text.  Variable to plot "BETA_STAR" or "BETA_HAT".  Only relevant if method includes 'Least_squares'
  heading = "Estimated family contributions to pooled samples" #Text. Title of thebar plot .  Only relevant if method includes 'Least_squares'
  plot.to.heading.height = 20 # Number. Height of the title relative to the height of thebar plot.  Only relevant if method includes 'Least_squares' 
  font.size.heading = 3       #Number. Font size ofbar plot heading.  Only relevant if method includes 'Least_squares'
  font.size.y.axis = 2        #Number. Font size ofbar plot y axis labels.  Only relevant if method includes 'Least_squares'
  font.size.x.axis = 2        #Number. Font size ofbar plot x axis labels.  Only relevant if method includes 'Least_squares'

    #Start checks############################################################################
  
  if(!skip.checks) {
    
    # method checks
    
    if(sum(method %in% c("Quantitative", "Discrete", "Exclusion", "Least_squares")) != length(method)) {
      stop("method must be a vector containing the words \'Quantitative\', \'Discrete\', \'Exclusion\', and or \'Least_squares\'")
    }  
    
    # discrete.method checks
    
    if(!(discrete.method %in% c("geno.probs", "assigned.genos"))) {
      stop("discrete.method must be \'geno.probs\' or \'assigned.genos\'")
    }
    
    # snp.dat.indiv checks 
    
    #Check that required headings are present in snp.dat.indiv  
    
    if(discrete.method  == "geno.probs") {
      
      if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B") %in% colnames(snp.dat.indiv)) != 4) {
        stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B")
      }  
      
      snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")]
      
      snp.dat.indiv$SAMPLE_ID   <- as.integer(snp.dat.indiv$SAMPLE_ID)
      snp.dat.indiv$SNP_ID      <- as.character(snp.dat.indiv$SNP_ID)  
      snp.dat.indiv$INTENSITY_A <- as.numeric(snp.dat.indiv$INTENSITY_A)  
      snp.dat.indiv$INTENSITY_B <- as.numeric(snp.dat.indiv$INTENSITY_B)  
      
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
      
    } 
    
    if(discrete.method  == "assigned.genos" & ("Quantitative" %in% method | "Least_squares" %in% method)) {
      
      if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "A_ALLELE", "B_ALLELE", "GENOTYPE") %in% colnames(snp.dat.indiv)) != 7) {
        stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B, A_ALLELE, B_ALLELE, GENOTYPE")
      }
      
      snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "A_ALLELE", "B_ALLELE", "GENOTYPE")] 
      
      snp.dat.indiv$SAMPLE_ID   <- as.integer(snp.dat.indiv$SAMPLE_ID)
      snp.dat.indiv$SNP_ID      <- as.character(snp.dat.indiv$SNP_ID)  
      snp.dat.indiv$INTENSITY_A <- as.numeric(snp.dat.indiv$INTENSITY_A)  
      snp.dat.indiv$INTENSITY_B <- as.numeric(snp.dat.indiv$INTENSITY_B) 
      snp.dat.indiv$A_ALLELE <- as.character(snp.dat.indiv$A_ALLELE)
      snp.dat.indiv$B_ALLELE <- as.character(snp.dat.indiv$B_ALLELE)
      snp.dat.indiv$GENOTYPE      <- as.character(snp.dat.indiv$GENOTYPE) 
      
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
    }
    
    if(discrete.method  == "assigned.genos" & !("Quantitative" %in% method | "Least_squares" %in% method)) {
      
      if(sum(c("SAMPLE_ID", "SNP_ID", "A_ALLELE", "B_ALLELE", "GENOTYPE") %in% colnames(snp.dat.indiv)) != 5) {
        stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, A_ALLELE, B_ALLELE, GENOTYPE")
      }
      
      snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "A_ALLELE", "B_ALLELE", "GENOTYPE")] 
      
      snp.dat.indiv$SAMPLE_ID   <- as.integer(snp.dat.indiv$SAMPLE_ID)
      snp.dat.indiv$SNP_ID      <- as.character(snp.dat.indiv$SNP_ID) 
      snp.dat.indiv$A_ALLELE <- as.character(snp.dat.indiv$A_ALLELE)
      snp.dat.indiv$B_ALLELE <- as.character(snp.dat.indiv$B_ALLELE)
      snp.dat.indiv$GENOTYPE <- as.character(snp.dat.indiv$GENOTYPE)   
      
      if(sum(!is.integer(snp.dat.indiv$SAMPLE_ID)) > 0) {
        stop("SAMPLE_ID in snp.dat.indiv must be an integer.  Also check for missing values.")
      }
      if(sum(!is.character(snp.dat.indiv$SNP_ID)) > 0) {
        stop("SNP_ID in snp.dat.indiv must be a character.  Check for missing values.")
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
    }
    
    if(sum(is.na(snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID")])) > 0) {
      stop("Check inputs in snp.dat.indiv  Are there missing SAMPLE_ID or SNP_ID?  Are columns of the correct class?")
    }
    
    #Check for duplicated records in snp.dat.pools
    indiv.snp <- paste(snp.dat.indiv$SAMPLE_ID,snp.dat.indiv$SNP_ID, sep=".")
    if(sum(duplicated(indiv.snp)) > 0) {
      stop("SAMPLE_ID and SNP_ID combinations are not unique in snp.dat.indiv  Delete duplicates or recode SAMPLE_ID.")
    }
    rm(indiv.snp)
    
    #order
    snp.dat.indiv <- snp.dat.indiv[order(snp.dat.indiv[,"SAMPLE_ID"]),]
    snp.dat.indiv <- snp.dat.indiv[order(snp.dat.indiv[,"SNP_ID"]),]
    
    # snp.dat.pools checks
    
    #Check that required headings are present in snp.dat.pools  
    
    if(discrete.method  == "geno.probs") {
      
      if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B") %in% colnames(snp.dat.pools)) != 4) {
        stop("snp.dat.pools input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B")
      }  
      
      snp.dat.pools <- snp.dat.pools[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")]
      
      snp.dat.pools$SAMPLE_ID   <- as.integer(snp.dat.pools$SAMPLE_ID)
      snp.dat.pools$SNP_ID      <- as.character(snp.dat.pools$SNP_ID)  
      snp.dat.pools$INTENSITY_A <- as.numeric(snp.dat.pools$INTENSITY_A)  
      snp.dat.pools$INTENSITY_B <- as.numeric(snp.dat.pools$INTENSITY_B)  
      
      if(sum(!is.integer(snp.dat.pools$SAMPLE_ID)) > 0) {
        stop("SAMPLE_ID in snp.dat.pools must be an integer.  Also check for missing values.")
      }
      if(sum(!is.character(snp.dat.pools$SNP_ID)) > 0) {
        stop("SNP_ID in snp.dat.pools must be a character.  Check for missing values.")
      }
      if(sum(!is.numeric(snp.dat.pools$INTENSITY_A)) > 0) {
        stop("INTENSITY_A in snp.dat.pools must be numeric.  Also check for missing values.")
      }
      if(sum(!is.numeric(snp.dat.pools$INTENSITY_B)) > 0) {
        stop("INTENSITY_B in snp.dat.pools must be numeric.  Also check for missing values.")
      }
      
    } 
    
    if(discrete.method  == "assigned.genos" & ("Quantitative" %in% method | "Least_squares" %in% method)) {
      
      if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "GENOTYPE") %in% colnames(snp.dat.pools)) != 5) {
        stop("snp.dat.pools input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B, GENOTYPE")
      }
      
      snp.dat.pools <- snp.dat.pools[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "GENOTYPE")] 
      
      snp.dat.pools$SAMPLE_ID   <- as.integer(snp.dat.pools$SAMPLE_ID)
      snp.dat.pools$SNP_ID      <- as.character(snp.dat.pools$SNP_ID)  
      snp.dat.pools$INTENSITY_A <- as.numeric(snp.dat.pools$INTENSITY_A)  
      snp.dat.pools$INTENSITY_B <- as.numeric(snp.dat.pools$INTENSITY_B)  
      snp.dat.pools$GENOTYPE      <- as.character(snp.dat.pools$GENOTYPE) 
      
      if(sum(!is.integer(snp.dat.pools$SAMPLE_ID)) > 0) {
        stop("SAMPLE_ID in snp.dat.pools must be an integer.  Also check for missing values.")
      }
      if(sum(!is.character(snp.dat.pools$SNP_ID)) > 0) {
        stop("SNP_ID in snp.dat.pools must be a character.  Check for missing values.")
      }
      if(sum(!is.numeric(snp.dat.pools$INTENSITY_A)) > 0) {
        stop("INTENSITY_A in snp.dat.pools must be numeric.  Also check for missing values.")
      }
      if(sum(!is.numeric(snp.dat.pools$INTENSITY_B)) > 0) {
        stop("INTENSITY_B in snp.dat.pools must be numeric.  Also check for missing values.")
      }
      if(sum(!is.character(snp.dat.pools$GENOTYPE)) > 0) {
        stop("GENOTYPE in snp.dat.pools must be a character.  Check for missing values.")
      }
    }
    
    if(discrete.method  == "assigned.genos" & !("Quantitative" %in% method | "Least_squares" %in% method)) {
      
      if(sum(c("SAMPLE_ID", "SNP_ID", "GENOTYPE") %in% colnames(snp.dat.pools)) != 3) {
        stop("snp.dat.pools input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, GENOTYPE")
      }
      
      snp.dat.pools <- snp.dat.pools[,c("SAMPLE_ID", "SNP_ID", "GENOTYPE")] 
      
      snp.dat.pools$SAMPLE_ID   <- as.integer(snp.dat.pools$SAMPLE_ID)
      snp.dat.pools$SNP_ID      <- as.character(snp.dat.pools$SNP_ID)  
      snp.dat.pools$GENOTYPE <- as.character(snp.dat.pools$GENOTYPE)   
      
      if(sum(!is.integer(snp.dat.pools$SAMPLE_ID)) > 0) {
        stop("SAMPLE_ID in snp.dat.pools must be an integer.  Also check for missing values.")
      }
      if(sum(!is.character(snp.dat.pools$SNP_ID)) > 0) {
        stop("SNP_ID in snp.dat.pools must be a character.  Check for missing values.")
      }
      if(sum(!is.character(snp.dat.pools$GENOTYPE)) > 0) {
        stop("GENOTYPE in snp.dat.pools must be a character.  Check for missing values.")
      }
    }
    
    if(sum(is.na(snp.dat.pools[,c("SAMPLE_ID", "SNP_ID")])) > 0) {
      stop("Check inputs in snp.dat.pools.  Are there missing SAMPLE_ID or SNP_ID?  Are columns of the correct class?")
    }
    
    #Check for duplicated records in snp.dat.pools
    indiv.snp <- paste(snp.dat.pools$SAMPLE_ID,snp.dat.pools$SNP_ID, sep=".")
    if(sum(duplicated(indiv.snp)) > 0) {
      stop("SAMPLE_ID and SNP_ID combinations are not unique in snp.dat.pools.  Delete duplicates or recode SAMPLE_ID.")
    }
    rm(indiv.snp)
    
    #order
    snp.dat.pools <- snp.dat.pools[order(snp.dat.pools[,"SAMPLE_ID"]),]
    snp.dat.pools <- snp.dat.pools[order(snp.dat.pools[,"SNP_ID"]),]
    
    # snp.param.indiv checks
    
    #Check that required headings are present in snp.param.indiv  
    if("Least_squares" %in% method &
       sum(c("WELCH_A", "WELCH_B") %in% colnames(snp.param.indiv)) != 2) {
      stop("If using the least squares method, the snp.param.indiv input must contain the following headings: WELCH_A and WELCH_B")
    }
    
    if((sum(c("Quantitative", "Discrete", "Exclusion") %in% method) > 0 &
        discrete.method == "geno.probs") |
        "Quantitative" %in% method) {
    if(sum(c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB", 
             "A_ALLELE", "B_ALLELE") %in% colnames(snp.param.indiv)) != 9) {  #"B_ALLELE_FREQ", 
      stop("snp.param.indiv input must contain the following headings: SNP_ID, MEAN_P_AA, SD_P_AA, MEAN_P_AB, SD_P_AB, MEAN_P_BB, SD_P_BB, A_ALLELE, and B_ALLELE") #B_ALLELE_FREQ,
    }
    
    snp.param.indiv$SNP_ID      <- as.character(snp.param.indiv$SNP_ID)  
    snp.param.indiv$MEAN_P_AA   <- as.numeric(snp.param.indiv$MEAN_P_AA)  
    snp.param.indiv$SD_P_AA     <- as.numeric(snp.param.indiv$SD_P_AA)  
    snp.param.indiv$MEAN_P_AB   <- as.numeric(snp.param.indiv$MEAN_P_AB)  
    snp.param.indiv$SD_P_AB     <- as.numeric(snp.param.indiv$SD_P_AB)  
    snp.param.indiv$MEAN_P_BB   <- as.numeric(snp.param.indiv$MEAN_P_BB)  
    snp.param.indiv$SD_P_BB     <- as.numeric(snp.param.indiv$SD_P_BB)  
    # snp.param.indiv$B_ALLELE_FREQ <- as.numeric(snp.param.indiv$B_ALLELE_FREQ)    
    snp.param.indiv$A_ALLELE    <- as.character(snp.param.indiv$A_ALLELE) 
    snp.param.indiv$B_ALLELE    <- as.character(snp.param.indiv$B_ALLELE) 
    
    if(sum(!is.character(snp.param.indiv$SNP_ID)) > 0) {
      stop("SNP_ID in snp.param.indiv must be a character.  Check for missing values.")
    }
    if(sum(!is.numeric(snp.param.indiv$MEAN_P_AA)) > 0) {
      stop("MEAN_P_AA in snp.param.indiv must be numeric.  Also check for missing values.")
    }
    if(sum(!is.numeric(snp.param.indiv$SD_P_AA)) > 0) {
      stop("SD_P_AA in snp.param.indiv must be numeric.  Also check for missing values.")
    }
    if(sum(!is.numeric(snp.param.indiv$MEAN_P_AB)) > 0) {
      stop("MEAN_P_AB in snp.param.indiv must be numeric.  Also check for missing values.")
    }
    if(sum(!is.numeric(snp.param.indiv$SD_P_AB)) > 0) {
      stop("SD_P_AB in snp.param.indiv must be numeric.  Also check for missing values.")
    }
    if(sum(!is.numeric(snp.param.indiv$MEAN_P_BB)) > 0) {
      stop("MEAN_P_BB in snp.param.indiv must be numeric.  Also check for missing values.")
    }
    if(sum(!is.numeric(snp.param.indiv$SD_P_BB)) > 0) {
      stop("SD_P_BB in snp.param.indiv must be numeric.  Also check for missing values.")
    }
    if(sum(!is.character(snp.param.indiv$A_ALLELE)) > 0) {
      stop("A_ALLELE in snp.param.indiv must be a character.  Check for missing values.")
    }
    if(sum(!is.character(snp.param.indiv$B_ALLELE)) > 0) {
      stop("B_ALLELE in snp.param.indiv must be a character.  Check for missing values.")
    }
    
    if(sum(is.na(snp.param.indiv)) > 0) {
      stop("Check inputs in snp.param.indiv  Are there missing values?  Are columns of the correct class?")
    }
    
    if("Least_squares" %in% method) {
      snp.param.indiv <- snp.param.indiv[,c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB", 
                                "WELCH_A", "WELCH_B",  "A_ALLELE", "B_ALLELE")]#"B_ALLELE_FREQ",
    } else {
      snp.param.indiv <- snp.param.indiv[,c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB", 
                                "A_ALLELE", "B_ALLELE")]#"B_ALLELE_FREQ", 
    }
    
    if(sum(snp.param.indiv[,c("SD_P_AA", "SD_P_AB", "SD_P_BB")] < 0) != 0) {
      stop("All values in the SD columns of snp.param.indiv must be greater than zero")
    }
    
    #Check that SNP are not duplicated
    if(sum(duplicated(snp.param.indiv[,"SNP_ID"])) > 0) {
      stop("SNP_IDs must not be duplicated in snp.param.indiv")
    }
    
    #Check that B_ALLELE_FREQ <=1 and >=0
    #  if(sum(snp.param.indiv[,"B_ALLELE_FREQ"] <=1 & snp.param.indiv[,"B_ALLELE_FREQ"] >= 0) != nrow(snp.param.indiv)) {
    #    stop("Data in the B_ALLELE_FREQ column must be between 0 and 1 inclusive")
    #  }       
    
    #Check that A_ALLELE and B_ALLELE data are one of A, C, G, or T
    if(sum(snp.param.indiv[,"A_ALLELE"] %in% c("A", "C", "G", "T") & 
           snp.param.indiv[,"B_ALLELE"] %in% c("A", "C", "G", "T")) !=  nrow(snp.param.indiv)) {
      stop("A_ALLELE and B_ALLELE data in snp.param.indiv must be \'A\', \'C\', \'G\', or \'T\'")
    }
    
    #Check that A_ALLELE is not the same as B_ALLELE
    if(sum(snp.param.indiv[,"A_ALLELE"] == snp.param.indiv[,"B_ALLELE"]) > 0) {
      stop("A_ALLELE and B_ALLELE must be different for each SNP in snp.param.indiv")
    }
    
    if("Least_squares" %in% method) {
      snp.param.indiv$WELCH_A <- as.numeric(snp.param.indiv$WELCH_A)  
      snp.param.indiv$WELCH_B <- as.numeric(snp.param.indiv$WELCH_B)  
      
      if(sum(!is.numeric(snp.param.indiv$WELCH_B)) > 0) {
        stop("WELCH_B in snp.param.indiv must be numeric.  Also check for missing values.")
      }
      if(sum(!is.numeric(snp.param.indiv$WELCH_B)) > 0) {
        stop("WELCH_B in snp.param.indiv must be numeric.  Also check for missing values.")
      }
    }
    
    #Check that there is only one nucleotide in column A_ALLELE for each SNP
    tmp.1 <-  unique(snp.param.indiv[,c("SNP_ID","A_ALLELE")])
    if(sum(unique(snp.param.indiv[,"SNP_ID"]) != tmp.1[,c("SNP_ID")]) > 0) {
      stop("There is more than one nucleotide in column A_ALLELE for at least one SNP in snp.param.indiv")
    }
    rm(tmp.1)
    
    #Check that there is only one nucleotide in column B_ALLELE for each SNP  
    tmp.1 <-  unique(snp.param.indiv[,c("SNP_ID","B_ALLELE")])
    if(sum(unique(snp.param.indiv[,"SNP_ID"]) != tmp.1[,c("SNP_ID")]) > 0) {
      stop("There is more than one nucleotide in column B_ALLELE for at least one SNP in snp.param.indiv")
    }
    rm(tmp.1)
    
    }
    
    #order
    snp.param.indiv <- snp.param.indiv[order(snp.param.indiv[,"SNP_ID"]),]
    
    # snp.param.pools checks
    
    if ((discrete.method == "geno.probs" | sum(c("Quantitative") %in% method) > 0) &
        is.null(snp.param.pools)) {
      stop("snp.param.pools must be specified in parent.assign.fun")
    }
    
    if ((sum(c("Quantitative", "Discrete", "Exclusion") %in% method) > 0 &
        discrete.method == "geno.probs") |
        "Quantitative" %in% method) {
      
      #Check that required headings are present in snp.param.pools  
      genotypes <- genotypes.fun(n.in.pools*2)
      
      if(sum(c("SNP_ID", paste("MEAN_P_", genotypes, sep=""), paste("SD_P_", genotypes, sep=""), "A_ALLELE", "B_ALLELE") %in% colnames(snp.param.pools)) != 
         (length(genotypes)*2 +3)) {
        stop(paste("snp.param.pools input (given the value of n.in.pools) must contain the following headings:", c("SNP_ID", 
                                                                                                                   paste("MEAN_P_", genotypes, sep=""), paste("SD_P_", genotypes, sep=""), "A_ALLELE", "B_ALLELE")))
      }
      
      snp.param.pools <- snp.param.pools[,c("SNP_ID", paste("MEAN_P_", genotypes, sep=""), paste("SD_P_", genotypes, sep=""), "A_ALLELE", "B_ALLELE")]
      
      snp.param.pools$SNP_ID   <- as.character(snp.param.pools$SNP_ID)  
      snp.param.pools$A_ALLELE <- as.character(snp.param.pools$A_ALLELE) 
      snp.param.pools$B_ALLELE <- as.character(snp.param.pools$B_ALLELE) 
      
      if(sum(!is.character(snp.param.pools$SNP_ID)) > 0) {
        stop("SNP_ID in snp.param.pools must be a character.  Check for missing values.")
      }
      if(sum(!is.character(snp.param.pools$A_ALLELE)) > 0) {
        stop("A_ALLELE in snp.param.pools must be a character.  Check for missing values.")
      }
      if(sum(!is.character(snp.param.pools$B_ALLELE)) > 0) {
        stop("B_ALLELE in snp.param.pools must be a character.  Check for missing values.")
      }
      
      for(genotype in genotypes) {
        tmp1 <- paste("MEAN_P_", genotype, sep="")
        tmp2 <- paste("SD_P_", genotype, sep="")
        snp.param.pools[,tmp1] <- as.numeric(snp.param.pools[,tmp1]) 
        snp.param.pools[,tmp2] <- as.numeric(snp.param.pools[,tmp2]) 
        rm(tmp1,tmp2)
      }
      
      if(sum(is.na(snp.param.pools)) > 0) {
        stop("Check inputs in snp.param.pools.  Are there missing values?  Are columns of the correct class?")
      }
      
      if(sum(snp.param.pools[,c(paste("SD_P_", genotypes, sep=""))] < 0) != 0) {
        stop("All values in the SD columns of snp.param.pools must be greater than zero")
      }
      
      #Check that SNP are not duplicated
      if(sum(duplicated(snp.param.pools[,"SNP_ID"])) > 0) {
        stop("SNP_IDs must no be duplicated in snp.param.pools")
      }
      
      #Check that A_ALLELE and B_ALLELE data are one of A, C, G, or T
      if(sum(snp.param.pools[,"A_ALLELE"] %in% c("A", "C", "G", "T") & 
             snp.param.pools[,"B_ALLELE"] %in% c("A", "C", "G", "T")) !=  nrow(snp.param.pools)) {
        stop("A_ALLELE and B_ALLELE data in snp.param.pools must be \'A\', \'C\', \'G\', or \'T\'")
      }
      
      #Check that A_ALLELE is not the same as B_ALLELE
      if(sum(snp.param.pools[,"A_ALLELE"] == snp.param.pools[,"B_ALLELE"]) > 0) {
        stop("A_ALLELE and B_ALLELE must be different for each SNP in snp.param.pools")
      }
      
      #Check SNP_ID and A_ALLELE/B_ALLELE combinations
      if(sum(unique(snp.param.pools[,c("SNP_ID", "A_ALLELE")]) == unique(snp.param.indiv[,c("SNP_ID", "A_ALLELE")])) 
         != (2*length(unique(snp.param.indiv$SNP_ID)))) {
        stop("SNP_ID and A_ALLELE combinations are different between snp.param.indiv and snp.param.pools")
      }
      
      #Check SNP_ID and A_ALLELE/B_ALLELE combinations
      if(sum(unique(snp.param.pools[,c("SNP_ID", "B_ALLELE")]) == unique(snp.param.indiv[,c("SNP_ID", "B_ALLELE")])) 
         != (2*length(unique(snp.param.indiv$SNP_ID)))) {
        stop("SNP_ID and B_ALLELE combinations are different between snp.param.indiv and snp.param.pools")
      }
      
      if(!nrow(unique(snp.param.pools[,c("SNP_ID", "A_ALLELE")])) == length(unique(snp.param.pools[,"SNP_ID"]))) {
        stop("Different A_ALLELEs listed for the same SNP_ID in snp.param.indiv and snp.param.pools")
      }
      
      if(!nrow(unique(snp.param.pools[,c("SNP_ID", "B_ALLELE")])) == length(unique(snp.param.pools[,"SNP_ID"]))) {
        stop("Different B_ALLELEs listed for the same SNP_ID in snp.param.indiv and snp.param.pools")
      }
      
    }
    
    #Check that there is only one nucleotide in column A_ALLELE for each SNP
    tmp.1 <-  unique(snp.param.pools[,c("SNP_ID","A_ALLELE")])
    if(sum(unique(snp.param.pools[,"SNP_ID"]) != tmp.1[,c("SNP_ID")]) > 0) {
      stop("There is more than one nucleotide in column A_ALLELE for at least one SNP in snp.param.pools")
    }
    rm(tmp.1)
    
    #Check that there is only one nucleotide in column B_ALLELE for each SNP  
    tmp.1 <-  unique(snp.param.pools[,c("SNP_ID","B_ALLELE")])
    if(sum(unique(snp.param.pools[,"SNP_ID"]) != tmp.1[,c("SNP_ID")]) > 0) {
      stop("There is more than one nucleotide in column B_ALLELE for at least one SNP in snp.param.pools")
    }
    rm(tmp.1)
    
    #order
    snp.param.pools <- snp.param.pools[order(snp.param.pools[,"SNP_ID"]),]
    
    #n.in.pools checks
    
    if(is.na(as.integer(n.in.pools))) {
      stop("n.in.pools must be an integer greater than zero")
    }
    
    n.in.pools <- as.integer(n.in.pools)
    
    if(n.in.pools < 1) {
      stop("n.in.pools must be an integer greater than zero")
    }
    
    # fams checks
    
    #Check that required headings are present in fams  
    if(sum(c("FAMILY_ID", "SIRE_ID", "DAM_ID") %in% colnames(fams)) != 3) {
      stop("fams input must be a data frame containing the following headings: FAMILY_ID, SIRE_ID, DAM_ID")
    }
    
    fams <- fams[,c("FAMILY_ID", "SIRE_ID", "DAM_ID")]
    
    fams$FAMILY_ID <- as.integer(fams$FAMILY_ID)
    fams$SIRE_ID   <- as.integer(fams$SIRE_ID)  
    fams$DAM_ID    <- as.integer(fams$DAM_ID)  
    
    if(sum(!is.integer(fams$FAMILY_ID)) > 0) {
      stop("FAMILY_ID in fams must be an integer.  Also check for missing values.")
    }
    if(sum(!is.integer(fams$SIRE_ID)) > 0) {
      stop("SIRE_ID in fams must be an integer.  Also check for missing values.")
    }
    if(sum(!is.integer(fams$DAM_ID)) > 0) {
      stop("DAM_ID in fams must be an integer.  Also check for missing values.")
    }
    
    if(sum(is.na(fams)) > 0) {
      stop("Check inputs in fams.  Are there missing values?  Are columns of the correct class?")
    }
    
    if(sum(fams < 1) > 0) {
      stop("FAMILY_ID, SIRE_ID and DAM_ID in fams must be integers and greater than 0")
    }
    
    if(!nrow(unique(fams[,c("SIRE_ID", "DAM_ID")])) == length(unique(fams[,c("FAMILY_ID")]))) {
      stop("There are multiple FAMILY_IDs with the same SIRE_ID and DAM_ID in fams")
    }
    
    if(sum(paste(fams[,"SIRE_ID"], fams[,"DAM_ID"], sep = "_") %in% 
           paste(fams[,"DAM_ID"], fams[,"SIRE_ID"], sep = "_")) > 0) {
      stop("There appear to be reciprocal families in fams.  These should be represented by only one row to avoid confusion.")
    }
    
    #Check that SIRE_ID and DAM_ID are present as SAMPLE_ID in snp.dat.indiv
    if(sum(!fams[,"SIRE_ID"] %in% snp.dat.indiv[,"SAMPLE_ID"]) > 0) {
      stop("All SIRE_IDs in fams must be present as SAMPLE_IDs in snp.dat.indiv")
    }
    
    #Check that SIRE_ID and DAM_ID are present as SAMPLE_ID in snp.dat.indiv
    if(sum(!fams[,"DAM_ID"] %in% snp.dat.indiv[,"SAMPLE_ID"]) > 0) {
      stop("All DAM_IDs in fams must be present as SAMPLE_IDs in snp.dat.indiv")
    }
    
    #order
    fams <- fams[order(fams[,"FAMILY_ID"]),]
    
    # fam.set.combns checks
    
    if((is.null(fam.set.combns.by.pool) & !is.null(fam.set.combns)) |
       (!is.null(fam.set.combns.by.pool) & is.null(fam.set.combns)) ) {
      stop("if either fam.set.combns.by.pool or fam.set.combns is NULL then both must be NULL")
    }
    
    if(!is.null(fam.set.combns)) {
    
    #Check that required headings are present in fam.set.combns  
    if(sum(c("FAM_SET_COMBN_ID", "FAM_SET_ID", "FAMILY_ID") %in% colnames(fam.set.combns)) != 3) {
      stop("fam.set.combns input must be a data frame containing the following headings: FAM_SET_COMBN_ID, FAM_SET_ID, FAMILY_ID")
    }
    
    fam.set.combns <- fam.set.combns[,c("FAM_SET_COMBN_ID", "FAM_SET_ID", "FAMILY_ID")]
    
    fam.set.combns$FAM_SET_COMBN_ID  <- as.integer(fam.set.combns$FAM_SET_COMBN_ID)
    fam.set.combns$FAM_SET_ID <- as.integer(fam.set.combns$FAM_SET_ID)  
    fam.set.combns$FAMILY_ID  <- as.integer(fam.set.combns$FAMILY_ID)  
    
    if(sum(!is.integer(fam.set.combns$FAM_SET_COMBN_ID)) > 0) {
      stop("FAM_SET_COMBN_ID in fam.set.combns must be an integer.  Also check for missing values.")
    }
    if(sum(!is.integer(fam.set.combns$FAM_SET_ID)) > 0) {
      stop("FAM_SET_ID in fam.set.combns must be an integer.  Also check for missing values.")
    }
    if(sum(!is.integer(fam.set.combns$FAMILY_ID)) > 0) {
      stop("FAMILY_ID in fam.set.combns must be an integer.  Also check for missing values.")
    }
    
    if(sum(is.na(fam.set.combns)) > 0) {
      stop("Check inputs in fam.set.combns.  Are there missing values?  Are columns of the correct class?")
    }
    
    if(sum(fam.set.combns[,c("FAM_SET_COMBN_ID", "FAM_SET_ID", "FAMILY_ID")] < 1) > 0) {
      stop("FAM_SET_COMBN_ID, FAM_SET_ID and FAMILY_ID in fam.set.combns must be integers and greater than 0")
    }
    
    tmp <- fam.set.combns
    tmp <- as.data.frame(table(tmp[,2:3]))
    tmp$FAM_SET_ID <- as.character(tmp$FAM_SET_ID)
    tmp$FAMILY_ID  <- as.character(tmp$FAMILY_ID)  
    tmp <- tmp[tmp[,"Freq"] != 0,]
    if(identical(!aggregate(tmp$Freq, by = list(tmp$FAM_SET_ID), na.rm=T, FUN = "mean")[2],
                 aggregate(tmp$Freq, by = list(tmp$FAM_SET_ID), na.rm=T, FUN = "max")[2])) {
      stop("FAM_SET_IDs present in multiple FAM_SET_COMBN_IDs do not allcontain the same FAMILY_IDs in fam.set.combns")
    }
    rm(tmp)
    
    #order
    fam.set.combns <- fam.set.combns[order(fam.set.combns[,"FAMILY_ID"]),]
    fam.set.combns <- fam.set.combns[order(fam.set.combns[,"FAM_SET_ID"]),]    
    fam.set.combns <- fam.set.combns[order(fam.set.combns[,"FAM_SET_COMBN_ID"]),]    
    }
    
    # fam.set.combns.by.pool checks
    
    if(!is.null(fam.set.combns.by.pool)) {
    
    #Check that required headings are present in fam.set.combns.by.pool  
    if(sum(c("SAMPLE_ID", "FAM_SET_COMBN_ID") %in% colnames(fam.set.combns.by.pool)) != 2) {
      stop("fam.set.combns.by.pool input must be a data frame containing the following headings: SAMPLE_ID, FAM_SET_COMBN_ID")
    }
    
    fam.set.combns.by.pool <- fam.set.combns.by.pool[,c("SAMPLE_ID", "FAM_SET_COMBN_ID")]
    
    fam.set.combns.by.pool$SAMPLE_ID <- as.integer(fam.set.combns.by.pool$SAMPLE_ID)      
    fam.set.combns.by.pool$FAM_SET_COMBN_ID  <- as.integer(fam.set.combns.by.pool$FAM_SET_COMBN_ID)
    
    if(sum(!is.integer(fam.set.combns.by.pool$SAMPLE_ID)) > 0) {
      stop("SAMPLE_ID in fam.set.combns.by.pool must be an integer.  Also check for missing values.")
    }
    
    if(sum(!is.integer(fam.set.combns.by.pool$FAM_SET_COMBN_ID)) > 0) {
      stop("FAM_SET_COMBN_ID in fam.set.combns.by.pool must be an integer.  Also check for missing values.")
    }
    
    if(sum(is.na(fam.set.combns.by.pool)) > 0) {
      stop("Check inputs in fam.set.combns.by.pool.  Are there missing values?  Are columns of the correct class?")
    }
    
    if(sum(fam.set.combns.by.pool[,c("SAMPLE_ID", "FAM_SET_COMBN_ID")] < 1) > 0) {
      stop("SAMPLE_ID and FAM_SET_COMBN_ID in fam.set.combns.by.pool must be integers and greater than 0")
    }
    
    fam.set.combns.by.pool <- fam.set.combns.by.pool[order(fam.set.combns.by.pool[,"SAMPLE_ID"]),]
    snp.dat.pools <- snp.dat.pools[order(snp.dat.pools[,"SAMPLE_ID"]),]
    if(!identical(as.character(unique(fam.set.combns.by.pool$SAMPLE_ID)), as.character(unique(snp.dat.pools$SAMPLE_ID)))) {
      stop("SAMPLE_ID in snp.dat.pools not the same as in fam.set.combns.by.pool")
    }
    
    fam.set.combns <- fam.set.combns[order(fam.set.combns[,"FAM_SET_COMBN_ID"]),]
    fam.set.combns.by.pool <- fam.set.combns.by.pool[order(fam.set.combns.by.pool[,"FAM_SET_COMBN_ID"]),]
    if(!identical(as.character(unique(fam.set.combns$FAM_SET_COMBN_ID)), as.character(unique(fam.set.combns.by.pool$FAM_SET_COMBN_ID)))) {
      stop("FAM_SET_COMBN_ID in fam.set.combns not the same as in fam.set.combns.by.pool")
    }
    
    #order
    fam.set.combns.by.pool <- fam.set.combns.by.pool[order(fam.set.combns.by.pool[,"SAMPLE_ID"]),]
    
    }
    
    # min.intensity checks
    
    if(is.na(as.numeric(min.intensity))) {
      stop("min.intensity must be a number")
    }
    
    min.intensity <- as.numeric(min.intensity)
    
    if(min.intensity < 0) {
      stop("min.intensity must be greater than or equal to zero")
    }
    
    # min.sd checks
    if(is.na(as.numeric(min.sd))) {
      stop("min.sd must be a number")
    }
    
    min.sd <- as.numeric(min.sd)
    
    if(min.sd < 0) {
      stop("min.sd must be greater than or equal to zero")
    }
    
    # snp.error.assumed checks
    
    if(!is.null(snp.error.assumed)) {
      
      if(is.data.frame(snp.error.assumed)) {
        #Check that required headings are present in fams  
        if(sum(c("SNP_ID", "SNP_ERROR_TILDE") %in% colnames(snp.error.assumed)) != 3) {
          stop("snp.error.assumed input must be a scalar or a data frame containing the following headings: SNP_ID, SNP_ERROR_TILDE")
        }
        
        snp.error.assumed <- snp.error.assumed[,c("SNP_ID", "SNP_ERROR_TILDE")]
        
        snp.error.assumed$SNP_ID <- as.character(snp.error.assumed$SNP_ID)
        snp.error.assumed$SNP_ERROR_TILDE   <- as.numeric(snp.error.assumed$SNP_ERROR_TILDE)  
        
        if(sum(is.na(snp.error.assumed)) > 0) {
          stop("Check inputs in snp.error.assumed  Are there missing values?  Are columns of the correct class?")
        }
        
        if(sum(snp.error.assumed$SNP_ERROR_TILDE < 0) > 0 |
           sum(snp.error.assumed$SNP_ERROR_TILDE > 1) > 0 ) {
          stop("SNP_ERROR_TILDE in snp.error.assumed must be a number between 0 and 1 inclusive")
        }
        
      } else {
        
        if(is.na(as.numeric(snp.error.assumed))) {
          stop("If snp.error.assumed is a scalar it must be between 0 and 1 inclusive")
        }
        
        snp.error.assumed <- as.numeric(snp.error.assumed)
        
        if(snp.error.assumed < 0 | snp.error.assumed > 1 ) {
          stop("If snp.error.assumed is a scalar it must be between 0 and 1 inclusive")
        }
        
      }
      
    }
    
    # snp.error.underlying checks
    
    if(sum(c("Quantitative", "Discrete", "Least_squares") %in% method) > 0) {
    
    if(is.null(snp.error.assumed)) {
      
      if(is.na(as.numeric(snp.error.underlying))) {
        stop("snp.error.underlying must be between 0 and 1 inclusive")
      }
      
      snp.error.underlying <- as.numeric(snp.error.underlying)
      
      if(snp.error.underlying. < 0 | snp.error.underlying > 1 ) {
        stop("snp.error.underlying must be between 0 and 1 inclusive")
      }
    }
    }
    
    # threshold.indiv checks
    
    if(sum(c("Discrete", "Exclusion") %in% method) > 0 &
       discrete.method == "geno.probs") {
      
      if(is.null(threshold.indiv)) {
        stop("threshold.indiv must be specified if method includes \'Discrete\' or \'Exclusion\' and discrete.method equals \'geno.probs\' ")
      }
      
      if(is.na(as.numeric(threshold.indiv))) {
        stop("threshold.indiv must be between 0 and 1 inclusive")
      }
      
      threshold.indiv <- as.numeric(threshold.indiv)
      
      if(threshold.indiv < 0 | threshold.indiv > 1 ) {
        stop("threshold.indiv must be between 0 and 1 inclusive")
      }
    }
    
    # threshold.pools checks
    
    if(sum(c("Discrete", "Exclusion") %in% method) > 0 &
       discrete.method == "geno.probs") {
      
      if(is.null(threshold.pools)) {
        stop("threshold.pools must be specified if method includes \'Discrete\' or \'Exclusion\'")
      }
      
      if(is.na(as.numeric(threshold.pools))) {
        stop("threshold.pools must be between 0 and 1 inclusive")
      }
      
      threshold.pools <- as.numeric(threshold.pools)
      
      if(threshold.pools < 0 | threshold.pools > 1 ) {
        stop("threshold.pools must be between 0 and 1 inclusive")
      }
    }
    
    
    if("Least_squares" %in% method) {
      
      # file.name checks 
      
      if(is.na(as.character(file.name))) {
        stop("file.name must be a string")
      }
      
      file.name <- as.character(file.name)
      
      # var checks   
      
      if(is.na(as.character(var))) {
        stop("var must be \'BETA_STAR\' or \'BETA_HAT\'")
      }
      
      var <- as.character(var)
      
      if(!var %in% c("BETA_STAR", "BETA_HAT")) {
        stop("var must be \'BETA_STAR\' or \'BETA_HAT\'")
      }
      
      # heading checks   
      
      if(is.na(as.character(heading))) {
        stop("heading must be a string")
      }
      
      heading <- as.character(heading)
      
      #plot.to.heading.height check 
      
      if(is.na(as.numeric(plot.to.heading.height))) {
        stop("plot.to.heading.height must be a number")
      }
      
      plot.to.heading.height <- as.numeric(plot.to.heading.height)
      
      if(plot.to.heading.height < 0) {
        stop("plot.to.heading.height must be greater than zero")
      }
      
      #font.size.heading check 
      
      if(is.na(as.numeric(font.size.heading))) {
        stop("font.size.heading must be a number")
      }
      
      font.size.heading <- as.numeric(font.size.heading)
      
      if(font.size.heading < 0) {
        stop("font.size.heading must be greater than zero")
      }
      
      #font.size.y.axis check 
      
      if(is.na(as.numeric(font.size.y.axis))) {
        stop("font.size.y.axis must be a number")
      }
      
      font.size.y.axis <- as.numeric(font.size.y.axis)
      
      if(font.size.y.axis < 0) {
        stop("font.size.y.axis must be greater than zero")
      }
      
      #font.size.y.axis check 
      
      if(is.na(as.numeric(font.size.x.axis))) {
        stop("font.size.x.axis must be a number")
      }
      
      font.size.x.axis <- as.numeric(font.size.x.axis)
      
      if(font.size.x.axis < 0) {
        stop("font.size.x.axis must be greater than zero")
      }
    }
    
    #Cross-data-frame checks
    
    #snp.dat.pools$SAMPLE_ID and snp.dat.indiv$SAMPLE_ID must all be different 
    if(sum(snp.dat.pools$SAMPLE_ID %in% snp.dat.indiv$SAMPLE_ID) > 0) {
      stop("SAMPLE_IDs in snp.dat.indiv must be different to SAMPLE_IDs in snp.dat.pools")
    }
    
    #    SNP_ID list must be the same in snp.dat.indiv, snp.dat.pools, snp.param.indiv, snp.param.pools (where relevant)
    tmp <- unique(snp.dat.indiv[order(snp.dat.indiv[,"SNP_ID"]),"SNP_ID"])
    if(!identical(tmp,unique(snp.dat.pools[order(snp.dat.pools[,"SNP_ID"]),"SNP_ID"]))) {
      stop("All SNP_IDs must be represented in both snp.dat.indiv and snp.dat.pools")
    }
    
    if ((sum(c("Discrete", "Exclusion") %in% method) > 0 &
         discrete.method == "geno.probs")|
        (sum(c("Quantitative", "Least_squares") %in% method)) > 0) {
          if(!identical(tmp,unique(snp.param.indiv$SNP_ID))) {
            stop("All SNP_IDs must be represented in each of snp.dat.indiv, snp.dat.pools, and snp.param.indiv")
          }
        }
        if ((sum(c("Discrete", "Exclusion") %in% method) > 0 &
             discrete.method == "geno.probs")|
            "Quantitative" %in% method) {
          if(!identical(tmp,unique(snp.param.pools$SNP_ID))) {
            stop("All SNP_IDs must be represented in each of snp.dat.indiv, snp.dat.pools, snp.param.indiv, snp.param.pools")
          }
        } 
    rm(tmp)
    
    #fams$SIRE_ID must be in snp.dat.indiv$SAMPLE_ID
    if(!sum(fams$SIRE_ID %in% snp.dat.indiv$SAMPLE_ID) == length(fams$SIRE_ID)) {
      stop ("Not all fams$SIRE_IDs are in snp.dat.indiv$SAMPLE_IDs")
    }
    
    #    fams$DAM_ID must be in snp.dat.indiv$SAMPLE_ID
    if(!sum(fams$DAM_ID %in% snp.dat.indiv$SAMPLE_ID) == length(fams$DAM_ID)) {
      stop ("Not all fams$DAM_IDs are in snp.dat.indiv$SAMPLE_IDs")
    }
    
    if(!is.null(fam.set.combns)){
    #fam.set.combns.by.pool$SAMPLE_ID must be in snp.dat.pools$SAMPLE_ID
    if(!sum(fam.set.combns.by.pool$SAMPLE_ID %in% snp.dat.pools$SAMPLE_ID) == length(fam.set.combns.by.pool$SAMPLE_ID)) {
      stop ("Not all fam.set.combns.by.pool$SAMPLE_IDs are in snp.dat.pools$SAMPLE_IDs")
    }
    
    #fam.set.combns$FAMILY_ID must be in fams$FAMILY_ID
    if(!sum(fam.set.combns$FAMILY_ID %in% fams$FAMILY_ID) == length(fam.set.combns$FAMILY_ID)) {
      stop ("Not all fam.set.combns$FAMILY_IDs are in fams$FAMILY_IDs")
    }
    }
    
    if(sum(c("Quantitative", "Discrete", "Least_squares") %in% method) > 0) {
      #snp.error.assumed and snp.error.underlying can't both be NULL
      if (is.null(snp.error.assumed) & is.null(snp.error.underlying)) {
        stop("Both snp.error.assumed and snp.error.underlying are NULL")
      }
    }
    
    #beta.min.ss
    if(!is.logical(beta.min.ss)) {
      stop("beta.min.ss must be either \'TRUE\' or \'FALSE\'")
    }
  }
  
  #End checks############################################################################
  
  # snp.dat.indiv <- snp.dat.indiv[,c("SNP_ID", "SAMPLE_ID", "GENOTYPE")]
  
  #Make outputs NULL
  phi.ij <- NULL
  Gij <- NULL
  flj.probs <- NULL
  flj.geno <- NULL
  parent.combns <- NULL
  lambda.kj <- NULL
  gkj <- NULL
  tclj.quant <- NULL
  tclj.discrete <- NULL
  tclj.ls  <- NULL             
  nlj.probs <- NULL
  nlj.geno <- NULL
  snp.error.probs <- NULL
  snp.error.geno <- NULL
  gklj.adj <- NULL
  dklj.adj <- NULL
  tclj.adj.quant <- NULL 
  tclj.adj.discrete <- NULL
  lod.duos.quant <- NULL
  lod.duos.discrete <- NULL
  logl.duos.quant <- NULL
  logl.duos.discrete <- NULL
  Dij  <- NULL                    
  dkj <- NULL                     
  mismatches  <- NULL             
  mismatches.by.snp <- NULL       
  # mismatches.by.snp.sample <- NULL
  # mismatch.snp.count <- NULL 
  Xl.mat <- NULL
  fkj.and.weight  <- NULL      
  beta     <- NULL            
  most.like.parents.quant <- NULL
  most.like.parents.discrete <- NULL
  most.like.parents.excl <- NULL 
  most.like.parents.excl.non.dup <- NULL 
  
  #if fam.set.combns and fam.set.combns.by.pool are NULL then generate data from 'fams' and  snp.dat.pools
  if(is.null(fam.set.combns) | is.null(fam.set.combns.by.pool)) {
    fam.set.combns <- data.frame(FAM_SET_COMBN_ID = 1,
      FAM_SET_ID = rep(1:n.in.pools, each = nrow(fams)),
      FAMILY_ID = rep(fams[,"FAMILY_ID"], n.in.pools))
    
    fam.set.combns.by.pool <- data.frame(SAMPLE_ID = unique(snp.dat.pools[,"SAMPLE_ID"]),
                                         FAM_SET_COMBN_ID = 1)
  }
  
  #if FAM_SET_COMBN_ID in fam.set.combns not in fam.set.combns.by.pool then remove
  fam.set.combns <- fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] %in% fam.set.combns.by.pool[,"FAM_SET_COMBN_ID"],] 

  #if FAMILY_ID in fams not in fam.set.combns then remove
  fams <- fams[fams[,"FAMILY_ID"] %in% fam.set.combns[,"FAMILY_ID"],] 
  
  #if individuals are not parents of families in fam.set.combns then delete from snp.dat.indiv and fams
  fams    <- fams[fams[,"FAMILY_ID"] %in% unique(fam.set.combns[,"FAMILY_ID"]),]
  snp.dat.indiv <- snp.dat.indiv[snp.dat.indiv[,"SAMPLE_ID"] %in% unique(c(fams[,"SIRE_ID"], fams[,"DAM_ID"])),]
  
  #run universal preliminary functions
  if(discrete.method == "geno.probs" | "Quantitative" %in% method | "Least_squares" %in% method ) { 
    phi.ij <- phi.ij.fun(snp.dat.indiv = snp.dat.indiv, 
                         snp.param.indiv = snp.param.indiv, 
                         fams = fams,
                         min.sd = min.sd,
                         min.intensity = min.intensity
    ) 
    
    Gij <- Gij.fun(phi.ij = phi.ij)
    
    flj.from.snp.dat <- NULL
    flj.probs    <- flj.from.parent.Gij.fun(Gij = Gij, fam.set.combns = fam.set.combns, fams = fams) 
  }  
  
  if(discrete.method == "assigned.genos") {   
    
  #  snp.dat.indiv$SNP_ID <- as.character(snp.dat.indiv$SNP_ID)
  #  snp.param.indiv$SNP_ID <- as.character(snp.param.indiv$SNP_ID)
  #  snp.dat.indiv <- left_join(snp.dat.indiv, unique(snp.param.indiv[,c("SNP_ID", "A_ALLELE", "B_ALLELE")]), by = "SNP_ID")
    
    flj.geno <- flj.from.snp.dat.fun(snp.dat.indiv = snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "A_ALLELE", "B_ALLELE", "GENOTYPE")], 
                                     fam.set.combns = fam.set.combns, 
                                     fams = fams)
  }  
  
  if (sum(c("Quantitative", "Discrete", "Exclusion") %in% method) > 0 ) {
    
    #get n.in.pools
  #  tmp <- colnames(snp.param.pools)[grep("MEAN_P_",colnames(snp.param.pools))][1]
  #  n.in.pools <- (nchar(tmp) - 7)/2
  #  rm(tmp)
    
    #get parent.combns
    tmp <- parent.combns.fun(fams = fams,
                             n.in.pools = n.in.pools, 
                             fam.set.combns = fam.set.combns)
    
    parent.combns         <- tmp$parent.combns
    parent.combns.by.fam.set.combn <- tmp$parent.combns.by.fam.set.combn
    
    rm(tmp)
    
  }       
  
  if((discrete.method == "geno.probs" & sum(c("Discrete", "Exclusion") %in% method) > 0 )|
     sum(c("Quantitative") %in% method) > 0 ) {          
    
    lambda.kj <- lambda.ij.fun(snp.dat.indiv = snp.dat.pools, 
                               snp.param.indiv = snp.param.pools, 
                               n.in.pools = n.in.pools,
                               min.sd = min.sd,
                               min.intensity = min.intensity
    )
    
    #get gkj
    rho.inv <- rho.inv.fun(n.in.pools)
    denominator <- rowSums(lambda.kj[,-c(1:2,ncol(lambda.kj))])
    
    gkj <- lambda.kj[,1:2]
    for(i in 1:nrow(rho.inv)) {
      gkj[,rho.inv[i,"GENOTYPE"]] <- rho.inv[i,"RHO_INV"] * lambda.kj[,2+i] / denominator
    }
    
    gkj <- gkj[order(gkj[,"SNP_ID"]),]
    gkj <- gkj[order(gkj[,"SAMPLE_ID"]),]
  }  
  
  if (sum(c("Quantitative", "Discrete", "Exclusion") %in% method) > 0 ) { 
    
    if(discrete.method == "geno.probs" | "Quantitative" %in% method) {  
      
      #compute nlj
      nlj.probs <- nlj.fun(flj = flj.probs,
                          n.in.pools = n.in.pools)  
      
      #get SNP error
      if(!is.null(snp.error.assumed)) {
        if(length(snp.error.assumed) == 1) {
          snp.error <- data.frame(SNP_ID = snp.param.indiv[,"SNP_ID"],
                                  SNP_ERROR_TILDE = snp.error.assumed)
        } else {
          snp.error <- snp.error.assumed
        }
      } else {
        snp.error.hat <- snp.error.fun(Gij = Gij,
                                       fams = fams, 
                                       parents.only = FALSE)
        snp.error <- snp.error.hat  
        snp.error[,"SNP_ERROR_TILDE"] <- snp.error[,"SNP_ERROR_HAT"]
        snp.error[snp.error[,"SNP_ERROR_TILDE"] < snp.error.underlying,"SNP_ERROR_TILDE"] <- snp.error.underlying
      }              
      snp.error.probs$SNP_ID <- as.character(snp.error.probs$SNP_ID)      
      snp.error.probs        <- snp.error[order(snp.error[,"SNP_ID"]),]
      
      rm(snp.error)
    }
    
    if("Discrete" %in% method & discrete.method == "assigned.genos") {   
      
      #compute nlj
      nlj.geno <- nlj.fun(flj = flj.geno,
                          n.in.pools = n.in.pools)  
      
      #get SNP error
      if(!is.null(snp.error.assumed)) {
        if(length(snp.error.assumed) == 1) {
          snp.error <- data.frame(SNP_ID = unique(snp.dat.indiv[,"SNP_ID"]),
                                  SNP_ERROR_TILDE = snp.error.assumed)
        } else {
          snp.error <- snp.error.assumed
        }
      } else {
        snp.error.hat <- snp.error.fun(Gij = Gij,
                                       fams = fams, 
                                       parents.only = FALSE)
        snp.error <- snp.error.hat  
        snp.error[,"SNP_ERROR_TILDE"] <- snp.error[,"SNP_ERROR_HAT"]
        snp.error[snp.error[,"SNP_ERROR_TILDE"] < snp.error.underlying,"SNP_ERROR_TILDE"] <- snp.error.underlying
      }              
      
      snp.error.geno <- snp.error[order(snp.error[,"SNP_ID"]),]
      rm(snp.error)
    }
  }
  
  if ("Quantitative" %in% method) {
    prelim.ml.quant <- prelim.ml.quant.fun(Gij = Gij,
                                           parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn,
                                           flj = flj.probs,
                                           gkj = gkj,
                                           snp.error = snp.error.probs,
                                           nlj = nlj.probs,
                                           fam.set.combns.by.pool =fam.set.combns.by.pool) 
    
    tclj.quant               = prelim.ml.quant$tclj
    miss.parent.count.quant = prelim.ml.quant$miss.parent.count
    gklj.adj     = prelim.ml.quant$gklj
    rm(prelim.ml.quant)
    
    ml <- ml.fun(g.d.klj.adj = gklj.adj,
                 tclj = tclj.quant,
                 snp.error = snp.error.probs,
                 nlj = nlj.probs,
                 miss.parent.count = miss.parent.count.quant,
                 parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn,
                 parent.combns = parent.combns,
                 meth = "quantitative")
    
    gklj.adj  <- ml$g.d.klj.adj
    tclj.adj.quant  <- ml$tclj.adj
    lod.duos.quant  <- ml$duos.lod
    logl.duos.quant <- ml$duos.logl
    most.like.parents.quant <- ml$most.like.parents
    rm(ml)
    
  }
  
  if (sum(c("Discrete", "Exclusion") %in% method) > 0 ) {
    
    if(discrete.method == "assigned.genos") {
      #     snp.dat.pools <- merge(snp.dat.pools, unique(snp.param.indiv[,c("SNP_ID", "A_ALLELE", "B_ALLELE")]), 
      #                             by = "SNP_ID", all.x = TRUE)
   #   snp.dat.pools$SNP_ID    <- as.character(snp.dat.pools$SNP_ID)
  #    snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  #    snp.dat.pools <- left_join(snp.dat.pools, unique(snp.param.indiv[,c("SNP_ID", "A_ALLELE", "B_ALLELE")]), by = "SNP_ID")
      
      prelim.ml.discrete <- prelim.ml.discrete.assigned.genos.fun(method = method,
                                                                  snp.dat.indiv = snp.dat.indiv,
                                                                  snp.dat.pools = snp.dat.pools,
                                                                  parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn,
                                                                  flj = flj.geno,
                                                                  snp.error = snp.error.geno,
                                                                  nlj = nlj.geno,
                                                                  fam.set.combns.by.pool =fam.set.combns.by.pool)
    } 
    
    if(discrete.method == "geno.probs") {
      prelim.ml.discrete <- prelim.ml.discrete.geno.probs.fun(
        Gij = Gij,
        threshold.indiv = threshold.indiv,
        snp.dat.indiv = snp.dat.indiv,
        parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn,
        flj = flj.probs,
        threshold.pools = threshold.pools,
        gkj = gkj,
        lambda.kj = lambda.kj,
        snp.error = snp.error.probs,
        nlj = nlj.probs,
        fam.set.combns.by.pool =fam.set.combns.by.pool)
    }
    
    Dij <- prelim.ml.discrete$Dij
    tclj.discrete <- prelim.ml.discrete$tclj
    miss.parent.count.discrete <- prelim.ml.discrete$miss.parent.count
    threshold.pools <- prelim.ml.discrete$threshold.pools
    dkj <- prelim.ml.discrete$dkj
    dklj.adj <- prelim.ml.discrete$g.d.klj.adj
    rm(prelim.ml.discrete)
  }  
  
  if ("Discrete" %in% method) {
    
    if(discrete.method == "geno.probs") {
      
      ml <- ml.fun(g.d.klj.adj = dklj.adj,
                   tclj = tclj.discrete,
                   snp.error = snp.error.probs,
                   nlj = nlj.probs,
                   miss.parent.count = miss.parent.count.discrete,
                   parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn,
                   parent.combns = parent.combns,
                   meth = "discrete")
      
      dklj.adj <- ml$g.d.klj.adj
      tclj.adj.discrete <- ml$tclj.adj
      lod.duos.discrete <- ml$duos.lod
      logl.duos.discrete <- ml$duos.logl
      most.like.parents.discrete <- ml$most.like.parents
      rm(ml)
      
    }
    
    if(discrete.method == "assigned.genos") {
      ml <- ml.fun(g.d.klj.adj = dklj.adj,
                   tclj = tclj.discrete,
                   snp.error = snp.error.geno,
                   nlj = nlj.geno,
                   miss.parent.count = miss.parent.count.discrete,
                   parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn,
                   parent.combns = parent.combns,
                   meth = "discrete")
      
      dklj.adj <- ml$g.d.klj.adj
      tclj.adj.discrete <- ml$tclj.adj
      lod.duos.discrete <- ml$duos.lod
      logl.duos.discrete <- ml$duos.logl
      most.like.parents.discrete <- ml$most.like.parents
      rm(ml)
    }
  }  
  
  if ("Exclusion" %in% method) {
    exclusion <- mismatches.fun(dkj = dkj,
                                tclj = tclj.discrete,
                                miss.parent.count = miss.parent.count.discrete,
                                parent.combns = parent.combns,
                                fam.set.combns.by.pool =fam.set.combns.by.pool) 
    
    mismatches               <- exclusion$duos.mismatches
    mismatches.by.snp        <- exclusion$duos.mismatches.by.snp
    mismatches.by.snp.sample <- exclusion$mismatches.by.snp.sample
    mismatch.snp.count       <- exclusion$mismatch.snp.count
    most.like.parents.excl   <- exclusion$most.like.parents
    most.like.parents.excl.non.dup   <- exclusion$most.like.parents.non.dup
    rm(exclusion)
  }  
  
  if ("Least_squares" %in% method) {
    
    #min.intensity.pools and n.by.pool functionality removed for general parent assign function  
    
    min.intensity.pools <- min.intensity 
    #   if(!beta.min.ss) {
    #      n.by.pool <- NULL
    #    } else {
    #      n.by.pool = data.frame(SAMPLE_ID = unique(fam.set.combns.by.pool$SAMPLE_ID),
    #                             N_INDIV = rep(n.in.pools, length(unique(fam.set.combns.by.pool$SAMPLE_ID))))
    #    }

    least.sq      <- ls.fun(fams = fams,
                            fam.set.combns = fam.set.combns,
                            fam.set.combns.by.pool =fam.set.combns.by.pool,
                            Gij = Gij,
                            flj = flj.probs,
                            snp.dat.pools = snp.dat.pools,
                            snp.param.indiv = snp.param.indiv,
                            min.intensity = min.intensity,
                            beta.min.ss = beta.min.ss)
    
    tclj.ls       <- least.sq$tclj.ls
    fkj.and.weight <- least.sq$fkj.and.weight
    Xl.mat         <- least.sq$Xl.mat
    beta          <- least.sq$beta
    rm(least.sq)
    
    if(!exists("running.sim")) { #don't plot if running simulations    
      wd <- getwd()
      dir.create(file.path(wd, "Results"), showWarnings = FALSE)
      setwd(file.path(wd, "Results"))
      
      for(samp in unique(beta$SAMPLE_ID)) {
        
        tmp.beta <- beta[beta[,"SAMPLE_ID"] == samp,]
        tmp.beta$FAMILY_ID <- as.character(tmp.beta$FAMILY_ID )
        
        bar.plot.fun(beta = tmp.beta,
                     file.name = samp,
                     var = var,
                     heading = heading,
                     plot.to.heading.height = plot.to.heading.height,
                     font.size.heading = font.size.heading,
                     font.size.y.axis = font.size.y.axis,
                     font.size.x.axis = font.size.x.axis)
      }
      setwd(wd)
    }
  }
  
  print(Sys.time())
  
  return(list(most.like.parents.quant = most.like.parents.quant,
              most.like.parents.discrete = most.like.parents.discrete,
              most.like.parents.excl = most.like.parents.excl,
              most.like.parents.excl.non.dup = most.like.parents.excl.non.dup,
              beta = beta,
              Dij = Dij,
              dkj = dkj,
              dklj.adj = dklj.adj,
              fkj.and.weight = fkj.and.weight,
              Gij = Gij,
              gkj = gkj,
              gklj.adj = gklj.adj,
              flj.probs = flj.probs,
              flj.geno = flj.geno,
              lambda.kj = lambda.kj,
              lod.duos.discrete = lod.duos.discrete,
              lod.duos.quant = lod.duos.quant,
              logl.duos.discrete = logl.duos.discrete,
              logl.duos.quant = logl.duos.quant,
              mismatches = mismatches,
              mismatches.by.snp = mismatches.by.snp,
              #    mismatches.by.snp.sample = mismatches.by.snp.sample,
              #    mismatch.snp.count = mismatch.snp.count,
              nlj.probs = nlj.probs,
              nlj.geno = nlj.geno,
              parent.combns = parent.combns,
              phi.ij = phi.ij,
              snp.error.probs = snp.error.probs,
              snp.error.geno = snp.error.geno,
              tclj.adj.quant = tclj.adj.quant,
              tclj.adj.discrete = tclj.adj.discrete,
              tclj.discrete = tclj.discrete,
              tclj.ls = tclj.ls,
              tclj.quant = tclj.quant,
              Xl.mat = Xl.mat
  ))
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

snp.param.pools.fun <- function(snp.param.indiv , #from snp.gen.param.fun
                                n.in.pools
) {
  
  #Estimates mean and standard deviation of allelic proportions for each possible genotype (combination of A and B alleles)
  #in DNA pools for each SNP 
  
  #Args##########################################
  
  # snp.param.indiv:   Data frame (output of snp.gen.param.fun) 
  #            1. SNP_ID    is the SNP identifier
  #            2. MEAN_P_AA   is the mean  of allelic proportion (homozygous allele A)
  #            3. SD_P_AA     is the standard deviation of allelic proporiton (homozygous allele A)
  #            4. MEAN_P_AB   is the mean  of allelic proportion (heterozygous)
  #            5. SD_P_AB     is the standard  deviation of allelic proporiton (heterozygous)
  #            6. MEAN_P_BB   is the mean of allelic proportion (homozygous allele B)
  #            7. SD_P_BB     is the standard deviation of allelic proporiton (homozygous allele B)
  #            8. A_ALLELE    is the base designated as the A allele
  #            9. B_ALLELE    is the base designated as the B allele
  
  #n.in.pools    Integer. Number of individuals in DNA pools. 
  
  #Returns##########################################
  
  # snp.param.pools  Data frame.  Number of columns depends on n.in.pools.  The following is for n.in.pools = 2.
  #            1. SNP_ID        is the SNP identifier
  #            2. MEAN_P_AAAA   is the mean of allelic proportion for individuals of genotype AAAA
  #            3. SD_P_AAAA     is the standard deviation of allelic proportion for individuals of genotype AAAA
  #            4. MEAN_P_AAAB  
  #            5. SD_P_AAAB 
  #            6. MEAN_P_AABB   
  #            7. SD_P_AABB 
  #            8. MEAN_P_ABBB   
  #            9. SD_P_ABBB 
  #            10. MEAN_P_BBBB   
  #            11. SD_P_BBBB 
  #            12. A_ALLELE    is the base designated as the A allele 
  #            13. B_ALLELE    is the base designated as the B allele
  
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

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

phi.ij.fun <- function(snp.dat.indiv, 
                       snp.param.indiv, 
                       fams,
                       min.sd = 0,
                       min.intensity = 0
) {  
  
  #Returns the elements ("PHI.AA", "PHI.AB", "PHI.BA", "PHI.BB") of the matrix shown in the 
  #top left of page 4 of Henshall et al. 2014.  Note that it is assumed that PHI.AB = PHI.BA.
  #Also returns the allelic proportion (pij) (see Equation 1 of Henshall et al. 2014). 
  
  #Retains parents only
  
  #Required functions:
  # pij.fun
  
  #Args##########################################
  # snp.dat.indiv: Data frame.  NOTE: only data for relevant offspring, sires and dams should be 
  #          included (i.e. data for non-relevant individuals and any pools should not be included 
  #          in this data frame).
  #            1. SAMPLE_ID  is the individual identifier
  #            2. SNP_ID    is the SNP identifier
  #            3. INTENSITY_A    is the area/intensity for allele A
  #            4. INTENSITY_B    is the area/intensity for allele B
  
  # snp.param.indiv:   Data frame (output of snp.gen.param.fun) 
  #            1. SNP_ID    is the SNP identifier
  #            2. MEAN_P_AA   is the mean  of allelic proportion (homozygous allele A)
  #            3. SD_P_AA     is the standard deviation of allelic proporiton (homozygous allele A)
  #            4. MEAN_P_AB   is the mean  of allelic proportion (heterozygous)
  #            5. SD_P_AB     is the standard  deviation of allelic proporiton (heterozygous)
  #            6. MEAN_P_BB   is the mean of allelic proportion (homozygous allele B)
  #            7. SD_P_BB     is the standard deviation of allelic proporiton (homozygous allele B)
  
  # min.sd       Number. Standard deviation of allelic proportion  fixed to this value if less 
  #              than it.
  
  # min.intensity      Number used in pij.fun. If sqrt((snp.dat.indiv$INTENSITY_A)^2 +
  #              (snp.dat.indiv$INTENSITY_B)^2) less than this value
  #              then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  #              Essentially removes observations that fall into the lower left of INTENSITY_A
  #              by INTENSITY_B scatter plot.
  
  #Returns##########################################
  # phi.ij:  Data frame 
  #            1.   SAMPLE_ID is the individual identifier
  #            2.   SNP_ID is the SNP identifier 
  #            3-6. AA_PHI, AB_PHI, BA_PHI, BB_PHI are the elements ("PHI.AA", "PHI.AB", "PHI.BA", 
  #                 "PHI.BB") of the matrix shown in the top left of page 4 of Henshall 
  #                 et al. 2014 
  #            7.   ALLELIC_PROP_INDIV is the allellic proportion (Equation 1 of Henshall 
  #                 et al. 2014) 
  
  print("Running phi.ij.fun")
  
  if(sum(c("FAMILY_ID", "SIRE_ID", "DAM_ID") %in% 
         colnames(fams)) != 3) {
    stop("fams input must be a data frame containing the following headings: FAMILY_ID, SIRE_ID, DAM_ID")
  }
  
  fams$FAMILY_ID    <- as.integer(fams$FAMILY_ID)
  fams$SIRE_ID   <- as.integer(fams$SIRE_ID)
  fams$DAM_ID    <- as.integer(fams$DAM_ID)
  fams <- fams[,c("FAMILY_ID", "SIRE_ID", "DAM_ID")]
  
  #If NA present as parents in fam then convert to 0
  fams[is.na(fams[,"SIRE_ID"]) ,"SIRE_ID"] <- "0"   
  fams[is.na(fams[,"DAM_ID"]) ,"DAM_ID"]   <- "0"  
  
  #identify parents
  parents <- unique(c(fams[,"SIRE_ID"],fams[,"DAM_ID"]))
  parents <- parents[parents != "0"]  
  
  #Check that all parents are in snp.dat.indiv.parents
  if(sum(!(parents %in% snp.dat.indiv[,"SAMPLE_ID"])) > 0) {
    stop("Not all parents in fams are present in snp.dat.indiv")
  }  
  
  #only retain parents in snp.dat.indiv
  snp.dat.indiv <- snp.dat.indiv[snp.dat.indiv[,"SAMPLE_ID"] %in% parents,]
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B") %in% 
         colnames(snp.dat.indiv)) != 4) {
    stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B")
  }
  
  if(sum(c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB") %in% 
         colnames(snp.param.indiv)) != 7) {
    stop("snp.param.indiv input must be a data frame containing the following headings: SNP_ID, MEAN_P_AA, SD_P_AA, MEAN_P_AB, SD_P_AB, MEAN_P_BB, SD_P_BB")
  }
  
  #Name columns and assign class
  # colnames(snp.dat.indiv) <- c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")
  snp.dat.indiv$SAMPLE_ID  <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$INTENSITY_A    <- as.numeric(snp.dat.indiv$INTENSITY_A)
  snp.dat.indiv$INTENSITY_B    <- as.numeric(snp.dat.indiv$INTENSITY_B)
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")]
  
  #colnames(snp.param.indiv) <- c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB")
  snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  snp.param.indiv$MEAN_P_AA   <- as.numeric(snp.param.indiv$MEAN_P_AA)
  snp.param.indiv$SD_P_AA     <- as.numeric(snp.param.indiv$SD_P_AA)
  snp.param.indiv$MEAN_P_AB   <- as.numeric(snp.param.indiv$MEAN_P_AB)
  snp.param.indiv$SD_P_AB     <- as.numeric(snp.param.indiv$SD_P_AB)
  snp.param.indiv$MEAN_P_BB   <- as.numeric(snp.param.indiv$MEAN_P_BB)
  snp.param.indiv$SD_P_BB     <- as.numeric(snp.param.indiv$SD_P_BB)
  snp.param.indiv <- snp.param.indiv[,c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB")]
  
  #Check for duplicated records in snp.dat.indiv
  indiv.snp <- paste(snp.dat.indiv$SAMPLE_ID,snp.dat.indiv$SNP_ID, sep=".")
  if(sum(duplicated(indiv.snp)) > 0) {
    stop("SAMPLE_ID and SNP_ID combinations are not unique in snp.dat.indiv.  Delete duplicates or recode SAMPLE_ID.")
  }
  rm(indiv.snp)
  
  # Check the list of SNPs the same in input files
  if(sum(unique(snp.dat.indiv[order(snp.dat.indiv[,"SNP_ID"]),"SNP_ID"]) != 
         unique(snp.param.indiv[order(snp.param.indiv[,"SNP_ID"]),"SNP_ID"]))>0) {
    stop("SNP identifiers do not match in snp.dat.indiv and snp.param.indiv")
  }
  
  #Get allelic proportion
  snp.allelic.prop <- pij.fun(snp.dat.indiv = snp.dat.indiv, min.intensity = min.intensity)
  print("Still running phi.ij.fun")
  
  #Rename columns
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "SAMPLE_ID"]    <- "SAMPLE_ID"
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "ALLELIC_PROP"] <- "ALLELIC_PROP_INDIV"  
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "INTENSITY"]    <- "INTENSITY_INDIV"  
  
  #merge data using cbind
  if(identical(snp.dat.indiv$SAMPLE_ID,snp.allelic.prop$SAMPLE_ID) &
     identical(snp.dat.indiv$SNP_ID,snp.allelic.prop$SNP_ID)) {
    snp.dat.indiv <- cbind(snp.dat.indiv, snp.allelic.prop[,!colnames(snp.allelic.prop) %in% c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")])
  } else {
    stop("SAMPLE_ID and SNP_ID columns of pij.fun output do not match those of snp.dat.indiv.  Not sure why.")
  }
  
  # Standard deviation of allelic proporiton fixed to min.sd if less than min.sd or is NA
  snp.param.indiv[snp.param.indiv[,"SD_P_AA"] < min.sd & !is.na(snp.param.indiv[,"SD_P_AA"]), "SD_P_AA"] <- min.sd
  snp.param.indiv[snp.param.indiv[,"SD_P_AB"] < min.sd & !is.na(snp.param.indiv[,"SD_P_AB"]), "SD_P_AB"] <- min.sd
  snp.param.indiv[snp.param.indiv[,"SD_P_BB"] < min.sd & !is.na(snp.param.indiv[,"SD_P_BB"]), "SD_P_BB"] <- min.sd
  
  snp.param.indiv[is.na(snp.param.indiv[,"SD_P_AA"]), "SD_P_AA"] <- min.sd
  snp.param.indiv[is.na(snp.param.indiv[,"SD_P_AB"]), "SD_P_AB"] <- min.sd
  snp.param.indiv[is.na(snp.param.indiv[,"SD_P_BB"]), "SD_P_BB"] <- min.sd
  
  #Ensure all SD are > 0 if mean is not NA
  if(sum(
    (
      !is.na(snp.param.indiv[,"MEAN_P_AA"]) &
      !is.na(snp.param.indiv[,"MEAN_P_AB"]) &
      !is.na(snp.param.indiv[,"MEAN_P_BB"])
    ) & 
    
    (
      (snp.param.indiv[,"SD_P_AA"] <= 0 & !is.na(snp.param.indiv[,"SD_P_AA"])) |
      (snp.param.indiv[,"SD_P_AB"] <= 0 & !is.na(snp.param.indiv[,"SD_P_AB"])) |
      (snp.param.indiv[,"SD_P_BB"] <= 0 & !is.na(snp.param.indiv[,"SD_P_BB"])) |
      is.na(snp.param.indiv[,"SD_P_AA"]) |
      is.na(snp.param.indiv[,"SD_P_AB"]) |
      is.na(snp.param.indiv[,"SD_P_BB"])
    )
  ) > 0) {
    stop("At least one element of SD_P_AA, SD_P_AB or SD_P_BB equals zero or is NA where P.AA, P.AB and P.BB are not NA in snp.param.indiv.  Must specify min.sd that is greater than 0 or modify snp.param.indiv.")
  }
  
  #Merge data
  # phi.ij <- merge(snp.dat.indiv, snp.param.indiv, by = "SNP_ID", all.x = TRUE)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  phi.ij <- left_join(snp.dat.indiv, snp.param.indiv, by = "SNP_ID")
  
  #Calculate G matrices - see start of page 4 of Henshall et al 2014
  phi.ij$AA_PHI <- dnorm(phi.ij$ALLELIC_PROP_INDIV,phi.ij$MEAN_P_AA,phi.ij$SD_P_AA)
  phi.ij$AB_PHI <- dnorm(phi.ij$ALLELIC_PROP_INDIV,phi.ij$MEAN_P_AB,phi.ij$SD_P_AB) / 2
  phi.ij$BA_PHI <- phi.ij$AB_PHI
  phi.ij$BB_PHI <- dnorm(phi.ij$ALLELIC_PROP_INDIV,phi.ij$MEAN_P_BB,phi.ij$SD_P_BB)
  
  #retain only relevant columns
  phi.ij <- phi.ij[,c("SAMPLE_ID", "SNP_ID", "AA_PHI", "AB_PHI", "BA_PHI", "BB_PHI", "ALLELIC_PROP_INDIV")]
  return(phi.ij)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

Gij.fun <- function(phi.ij, round.Gij = NULL) {
  
  #Returns########################################## the elements of the Gij matrix ("AA_GENO_PROB", " AB_GENO_PROB", " BA_GENO_PROB", "BB_GENO_PROB") 
  #and Tij vector ("A_TRANS_PROB", "B_TRANS_PROB") shown in the top left of page 4 
  #of Henshall et al. 2014.  
  
  #Args##########################################
  # phi.ij:  Data frame (output of phi.ij.fun)
  #            1.   SAMPLE_ID is the individual identifier
  #            2.   SNP_ID is the SNP identifier 
  #            3-6. AA_PHI, AB_PHI, BA_PHI, BB_PHI are the elements ("PHI.AA", "PHI.AB", "PHI.BA", 
  #                 "PHI.BB") of the matrix shown in the top left of page 4 of Henshall 
  #                 et al. 2014 
  #            7.   ALLELIC_PROP_INDIV is the allellic proportion (Equation 1 of Henshall 
  #                 et al. 2014) 
  
  #round.Gij: Integer.  Generally NULL. Rounds GENO_PROB values to specified number of decimal points.
  #           NOTE: when the round.Gij is specified (i.e. not NULL) there is a high probability that Equation 3 will equal 0
  #           (which essentially suggests a mismatch) and will result in infinite likelihood ratios.
  
  #Returns##########################################
  # Gij:     Data frame 
  #              1.   SAMPLE_ID is the individual identifier 
  #              2.   SNP_ID is the SNP identifier
  #              3-6. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the Gij matrix (see the 
  #                   top left of page 4 of Henshall et al. 2014
  #              7-8. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #                   for alleles A and B respectively computed from Gij (i.e. the elements of 
  #                   the transmission (Tij) vector, Equation 2 of Henshall et al. 2014).
  #              9.   ALLELIC_PROP_INDIV is the allellic proportion (Equation 1 of Henshall 
  #                   et al. 2014).  Retained from the input data frame.
  
  print("Running Gij.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "AA_PHI", "AB_PHI", "BA_PHI", "BB_PHI", "ALLELIC_PROP_INDIV") %in% 
         colnames(phi.ij)) != 7) {
    stop("phi.ij input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, AA_PHI, AB_PHI, BA_PHI, BB_PHI, ALLELIC_PROP_INDIV")
  }
  
  #Change column names and checks
  phi.ij$SAMPLE_ID  <- as.integer(phi.ij$SAMPLE_ID)
  phi.ij$SNP_ID    <- as.character(phi.ij$SNP_ID)
  phi.ij$AA_PHI <- as.numeric(phi.ij$AA_PHI)
  phi.ij$AB_PHI <- as.numeric(phi.ij$AB_PHI)
  phi.ij$BA_PHI <- as.numeric(phi.ij$BA_PHI)
  phi.ij$BB_PHI <- as.numeric(phi.ij$BB_PHI)
  phi.ij$ALLELIC_PROP_INDIV <- as.numeric(phi.ij$ALLELIC_PROP_INDIV)
  phi.ij <- phi.ij[,c("SAMPLE_ID", "SNP_ID", "AA_PHI", "AB_PHI", "BA_PHI", "BB_PHI", "ALLELIC_PROP_INDIV")]
  
  Gij <- phi.ij
  
  #G matrix
  phi.sum         <- Gij$AA_PHI + Gij$AB_PHI + Gij$BA_PHI + Gij$BB_PHI
  Gij$AA_GENO_PROB <- Gij$AA_PHI / phi.sum 
  Gij$AB_GENO_PROB <- Gij$AB_PHI / phi.sum
  Gij$BA_GENO_PROB <- Gij$BA_PHI / phi.sum
  Gij$BB_GENO_PROB <- Gij$BB_PHI / phi.sum 
  rm(phi.sum)
  
  
  if(!is.null(round.Gij)) {  
    
    #Must sum to 1
    #https://www.r-bloggers.com/round-values-while-preserve-their-rounded-sum-in-r/
    round_preserve_sum <- function(x) {
      
      x_orig <- x
      x[is.na(x)] <- 0.25
      
      up <- 10 ^ round.Gij
      x <- x * up
      y <- floor(x)
      indices <- tail(order(x-y), round(sum(x)) - sum(y))
      if(length(indices) > 0) {
        y[indices] <- y[indices] + 1
      }
      z <- y / up
      
      z[is.na(x_orig)] <- NA
      z
    }
    
    #Apply function by row
    Gij[,c("AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")] <-  
      apply(as.matrix(Gij[,c("AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")]), 2, round_preserve_sum)
    
    #Get some odd results (i.e. don't sum exactly to one)
    G.sum         <- Gij$AA_GENO_PROB + Gij$AB_GENO_PROB + Gij$BA_GENO_PROB + Gij$BB_GENO_PROB 
    Gij$AA_GENO_PROB <- Gij$AA_GENO_PROB / G.sum 
    Gij$AB_GENO_PROB <- Gij$AB_GENO_PROB / G.sum
    Gij$BA_GENO_PROB <- Gij$BA_GENO_PROB / G.sum 
    Gij$BB_GENO_PROB <- Gij$BB_GENO_PROB / G.sum
    rm(G.sum)
    
    #Rerun - Apply function by row
    Gij[,c("AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")] <-  
      apply(as.matrix(Gij[,c("AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")]), 2, round_preserve_sum)
    
    #Get some odd results (i.e. don't sum exactly to one)
    G.sum         <- Gij$AA_GENO_PROB + Gij$AB_GENO_PROB + Gij$BA_GENO_PROB + Gij$BB_GENO_PROB 
    Gij$AA_GENO_PROB <- Gij$AA_GENO_PROB / G.sum 
    Gij$AB_GENO_PROB <- Gij$AB_GENO_PROB / G.sum
    Gij$BA_GENO_PROB <- Gij$BA_GENO_PROB / G.sum 
    Gij$BB_GENO_PROB <- Gij$BB_GENO_PROB / G.sum
    rm(G.sum)
    
  }
  
  #Equation 2 of Henshall et al. 2014 - T vector
  Gij$A_TRANS_PROB <- (Gij$AA_GENO_PROB * 2 + Gij$AB_GENO_PROB + Gij$BA_GENO_PROB)/2  #(sum row 1 + sum col 1) / 2
  Gij$B_TRANS_PROB <- (Gij$BB_GENO_PROB * 2 + Gij$AB_GENO_PROB + Gij$BA_GENO_PROB)/2  #(sum row 2 + sum col 2) / 2
  
  #Get some odd results (i.e. don't sum exactly to one)
  G.sum         <- Gij$A_TRANS_PROB + Gij$B_TRANS_PROB 
  Gij$A_TRANS_PROB <- Gij$A_TRANS_PROB / G.sum 
  Gij$B_TRANS_PROB <- Gij$B_TRANS_PROB / G.sum
  rm(G.sum)
  
  #retain only relevant columns
  Gij <- Gij[,c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB",
                "A_TRANS_PROB", "B_TRANS_PROB", "ALLELIC_PROP_INDIV")]
  return(Gij)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

flj.from.parent.Gij.fun <- function(Gij, fam.set.combns, fams) {
  
  #Generates allele frequencies by SNP from means of G matrices of parents. Used in early drafts of Henshall et al. 2014.
  
  #Args##########################################
  
  #Gij: Data frame.  Output of Gij.fun 
  #          1.   SAMPLE_ID is the individual identifier (assumed to be offspring if not 
  #               listed as a sire or dam in 'fams' input)
  #          2.   SNP_ID is the SNP identifier
  #          3-6. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the Gij  matrix (see the 
  #               top left of page 4 of Henshall et al. 2014
  #          7-8. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #               for alleles A and B respectively computed from Gij (i.e. the elements of 
  #               the transmission (Tij) vector, Equation 2 of Henshall et al. 2014).
  
  # fam.set.combns:    Data frame.  NOTE: to define half sib families, the unknown parent must be NA or 0
  #                       1. FAM_SET_COMBN_ID (integer) itentifer of unique combination of FAM_SET_ID
  #                       2. FAM_SET_ID (integer) is the family aggregate identifier (integer > 0) 
  #                       3. FAMILY_ID (integer) is the family identifier  
  
  # fams:    Data frame.  NOTE: to define half sib families, the unknown parent must be NA or 0
  #          1. FAMILY_ID is the family identifier
  #          2. SIRE_ID is the sire identifier 
  #          3. DAM_ID is the dam identifier 
  
  #Returns##########################################
  
  # flj: Data frame.  Output of flj.from.parent.Gij.fun 
  #          1.   SNP_ID is the snp identifier
  #          2.     FAM_SET_COMBN_ID
  #          . AA_GENO_PROB_MEAN, AB_GENO_PROB_MEAN, BA_GENO_PROB_MEAN, BB_GENO_PROB_MEAN: means of parental
  #               G matrices
  #          . A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #               for alleles A and B respectively (i.e. the elements of the flj vector) computed from average parental G matrices 
  #               (i.e. from AA_GENO_PROB_MEAN, AB_GENO_PROB_MEAN, BA_GENO_PROB_MEAN, BB_GENO_PROB_MEAN).  Estimates of A and B allele 
  #               frequencies in the parental 'population'.
  #          . AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the flj matrix (see the 
  #               bottom left of page 4 of Henshall et al. 2014
  #          . SAMPLE_ID is the individual identifier (all values are zero and represent unknown individuals)
  
  print("Running flj.from.parent.Gij.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB") %in% 
         colnames(Gij)) != 4) {
    stop("Gij input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, A_TRANS_PROB, B_TRANS_PROB")
  }
  
  if(sum(c("FAM_SET_COMBN_ID", "FAM_SET_ID", "FAMILY_ID") %in% 
         colnames(fam.set.combns)) != 3) {
    stop("fam.set.combns input must be a data frame containing the following headings: FAM_SET_COMBN_ID FAM_SET_ID FAMILY_ID")
  }
  
  if(sum(c("FAMILY_ID", "SIRE_ID", "DAM_ID") %in% 
         colnames(fams)) != 3) {
    stop("fams input must be a data frame containing the following headings: FAMILY_ID, SIRE_ID, DAM_ID")
  }
  
  #Name columns and assign class
  #colnames(Gij)     <- c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB", 
  #                                 "A_TRANS_PROB", "B_TRANS_PROB")
  Gij$SAMPLE_ID     <- as.integer(Gij$SAMPLE_ID)
  Gij$SNP_ID       <- as.character(Gij$SNP_ID)
  Gij$A_TRANS_PROB <- as.numeric(Gij$A_TRANS_PROB)
  Gij$B_TRANS_PROB <- as.numeric(Gij$B_TRANS_PROB)
  Gij <- Gij[,c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")]
  
  fam.set.combns$FAM_SET_COMBN_ID <- as.integer(fam.set.combns$FAM_SET_COMBN_ID)
  fam.set.combns$FAM_SET_ID       <- as.integer(fam.set.combns$FAM_SET_ID)
  fam.set.combns$FAMILY_ID        <- as.integer(fam.set.combns$FAMILY_ID)
  
  #colnames(fams) <- c("FAMILY_ID", "SIRE_ID", "DAM_ID")
  fams$FAMILY_ID    <- as.integer(fams$FAMILY_ID)
  fams$SIRE_ID   <- as.integer(fams$SIRE_ID)
  fams$DAM_ID    <- as.integer(fams$DAM_ID)
  fams <- fams[,c("FAMILY_ID", "SIRE_ID", "DAM_ID")]
  
  #Check that probabilities between 0 and 1
  if(
    sum(
      Gij$A_TRANS_PROB > 1 |
      Gij$B_TRANS_PROB > 1 |
      
      Gij$A_TRANS_PROB < 0 |
      Gij$B_TRANS_PROB < 0 
      , na.rm = T) != 0
  ) {
    stop("Probabilities in Gij must be between 0 and 1 inclusive")
  }
  
  #Check that "A_TRANS_PROB", "B_TRANS_PROB" sum to 1 in Gij
  if(
    sum(
      round((Gij$A_TRANS_PROB + Gij$B_TRANS_PROB),5) != 1.0, na.rm = T
    ) != 0
  ) {
    stop("A_TRANS_PROB + B_TRANS_PROB must equal 1 in all rows of Gij")
  }
  
  #Ensure that SAMPLE_ID column contains no 0s or NA
  if(
    sum(
      (is.na(Gij$SAMPLE_ID) |
       Gij$SAMPLE_ID == 0), na.rm = T
    ) != 0
  ) {
    stop("SAMPLE_ID in Gij cannot be NA or 0")
  } 
  
  #If NA present as parents in fam then convert to 0
  fams[is.na(fams[,"SIRE_ID"]) ,"SIRE_ID"] <- 0   
  fams[is.na(fams[,"DAM_ID"]) ,"DAM_ID"]   <- 0  
  
  flj.from.parent.G <- NULL
  
  #cycle through FAM_SET_COMBN_ID
  
  for(fam.set.combn in unique(fam.set.combns$FAM_SET_COMBN_ID)) {
    
    #identify families in fam.set.combn
    
    tmp.fams <- fams[fams[,"FAMILY_ID"] %in% fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn,"FAMILY_ID"],]
    
    #get parent subset of Gij
    
    #identify parents
    parents <- data.frame(SAMPLE_ID = c(tmp.fams[,"SIRE_ID"],tmp.fams[,"DAM_ID"]))
    parents$SAMPLE_ID <- parents[parents$SAMPLE_ID != 0,]  
    parents$SAMPLE_ID <- as.integer(as.character(parents$SAMPLE_ID))
    
    #separate parents and offspring in Gij
    G.T.parents <- Gij[Gij[,"SAMPLE_ID"] %in% parents,]
    G.T.parents <- left_join(parents, Gij, by = "SAMPLE_ID")
    
    #get vector of allele probabilities of each SNP estimated from all parents
    a.mean.parents <- aggregate(G.T.parents$A_TRANS_PROB, by = list(G.T.parents$SNP_ID), 
                                na.rm=T, FUN = "mean")   
    colnames(a.mean.parents) <- c("SNP_ID", "A_TRANS_PROB")
    
    b.mean.parents <- aggregate(G.T.parents$B_TRANS_PROB, by = list(G.T.parents$SNP_ID), 
                                na.rm=T, FUN = "mean")   
    colnames(b.mean.parents) <- c("SNP_ID", "B_TRANS_PROB")
    
    #  tmp.flj.from.parent.G <- merge(a.mean.parents, b.mean.parents, by = "SNP_ID", all = TRUE)
    a.mean.parents$SNP_ID    <- as.character(a.mean.parents$SNP_ID)
    b.mean.parents$SNP_ID    <- as.character(b.mean.parents$SNP_ID)
    tmp.flj.from.parent.G <- inner_join(a.mean.parents, b.mean.parents, by = "SNP_ID")
    
    #Add rows for missing parents to G.T.parents  (see bottom left page 4 of Henshall et al. 2014)
    tmp.flj.from.parent.G$SAMPLE_ID <- as.integer(0)
    
    tmp.flj.from.parent.G$FAM_SET_COMBN_ID <- fam.set.combn
    
    flj.from.parent.G <- rbind(flj.from.parent.G,tmp.flj.from.parent.G)
  }
  
  rm(tmp.fams, tmp.flj.from.parent.G, fam.set.combn)
  flj.from.parent.G <- flj.from.parent.G[,c("FAM_SET_COMBN_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB", "SAMPLE_ID")]
  return(flj.from.parent.G)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

flj.from.snp.dat.fun <- function(snp.dat.indiv, fam.set.combns, fams) {
  
  #Generates SNP allele frequencies from genotype calls for parents only
  
  #Required functions:
  # pij.fun, snp.gen.param.fun
  
  #Args##########################################
  # snp.dat.indiv: Data frame
  #              SAMPLE_ID is the individual identifier
  #              SNP_ID   is the SNP identifier
  #              A_ALLELE  is the base represented by allele A
  #              B_ALLELE  is the base represented by allele B
  #              GENOTYPE is the SNP genotype call
  
  # fam.set.combns:    Data frame.  NOTE: to define half sib families, the unknown parent must be NA or 0
  #                       1. FAM_SET_COMBN_ID (integer) itentifer of unique combination of FAM_SET_ID
  #                       2. FAM_SET_ID (integer) is the family aggregate identifier (integer > 0) 
  #                       3. FAMILY_ID (integer) is the family identifier  
  
  # fams: Data frame.  NOTE: to define half sib families, the unknown parent must be NA or 0
  #          1. FAMILY_ID is the family identifier
  #          2. SIRE_ID is the sire identifier 
  #          3. DAM_ID is the dam identifier 
  
  #Returns##########################################
  #   Data frame. See "Estimation of SNP sepecific parameters" page 3 of Henshall et al 2014 
  #              SNP_ID        is the SNP identifier, 
  #              N_AA          is the count homozygotes for allele A
  #              N_AB          is the count heterozygotes
  #              N_BB          is the count homozygotes for allele B
  #              A_ALLELE_FREQ is the A allele frequency derived from counts
  #              B_ALLELE_FREQ is the B allele frequency derived from counts
  
  print("Running flj.from.snp.dat.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID",  "A_ALLELE", "B_ALLELE", "GENOTYPE") %in% colnames(snp.dat.indiv)) != 5) {
    stop("snp.dat.indiv input for flj.from.snp.dat.fun must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, A_ALLELE, B_ALLELE, GENOTYPE")
  }
  
  if(sum(c("FAM_SET_COMBN_ID", "FAM_SET_ID", "FAMILY_ID") %in% 
         colnames(fam.set.combns)) != 3) {
    stop("fam.set.combns input must be a data frame containing the following headings: FAM_SET_COMBN_ID FAM_SET_ID FAMILY_ID")
  }
  
  if(sum(c("FAMILY_ID", "SIRE_ID", "DAM_ID") %in% 
         colnames(fams)) != 3) {
    stop("fams input must be a data frame containing the following headings: FAMILY_ID, SIRE_ID, DAM_ID")
  }
  
  #Retain columns and assign class
  #colnames(snp.dat.indiv) <- c("SAMPLE_ID", "SNP_ID", "A_ALLELE", "B_ALLELE", "GENOTYPE")
  snp.dat.indiv$SAMPLE_ID  <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$A_ALLELE   <- as.character(snp.dat.indiv$A_ALLELE)
  snp.dat.indiv$B_ALLELE   <- as.character(snp.dat.indiv$B_ALLELE)
  snp.dat.indiv$GENOTYPE  <- as.character(snp.dat.indiv$GENOTYPE)
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "A_ALLELE", "B_ALLELE", "GENOTYPE")]
  
  fam.set.combns$FAM_SET_COMBN_ID <- as.integer(fam.set.combns$FAM_SET_COMBN_ID)
  fam.set.combns$FAM_SET_ID       <- as.integer(fam.set.combns$FAM_SET_ID)
  fam.set.combns$FAMILY_ID        <- as.integer(fam.set.combns$FAMILY_ID)
  
  #colnames(fams) <- c("FAMILY_ID", "SIRE_ID", "DAM_ID")
  fams$FAMILY_ID    <- as.integer(fams$FAMILY_ID)
  fams$SIRE_ID   <- as.integer(fams$SIRE_ID)
  fams$DAM_ID    <- as.integer(fams$DAM_ID)
  fams <- fams[,c("FAMILY_ID", "SIRE_ID", "DAM_ID")]
  
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
  
  #If NA present as parents in fam then convert to 0
  fams[is.na(fams[,"SIRE_ID"]) ,"SIRE_ID"] <- "0"   
  fams[is.na(fams[,"DAM_ID"]) ,"DAM_ID"]   <- "0"  
  
  #Dummy data required to run snp.gen.param.fun
  snp.dat.indiv$INTENSITY_A    <- 10000*runif(nrow(snp.dat.indiv)) #random numbers
  snp.dat.indiv$INTENSITY_B    <- 10000*runif(nrow(snp.dat.indiv))  
  
  flj.from.snp.dat <- NULL
  
  for(fam.set.combn in unique(fam.set.combns$FAM_SET_COMBN_ID)) {
    #identify families in fam.set.combn
    tmp.fams <- fams[fams[,"FAMILY_ID"] %in% fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn,"FAMILY_ID"],]
    
    #identify parents
    parents <- data.frame(SAMPLE_ID = c(tmp.fams[,"SIRE_ID"],tmp.fams[,"DAM_ID"]))
    parents$SAMPLE_ID <- parents[parents$SAMPLE_ID != 0,]  
    parents$SAMPLE_ID <- as.integer(as.character(parents$SAMPLE_ID))
    #generate dummy SAMPLE_ID to allow snp.gen.param.fun to run (doesn't allow duplicate SAMPLE_ID/SNP_ID combinations)
    parents$DUMMY_SAMPLE_ID <- 1:nrow(parents) 
    #Get allele frequencies from counts using snp.gen.param.fun
    
    #duplciate snp data if a parent represented more than once
    tmp.snp.dat.indiv <- left_join(parents, snp.dat.indiv, by = "SAMPLE_ID")
    tmp.snp.dat.indiv$SAMPLE_ID <- tmp.snp.dat.indiv$DUMMY_SAMPLE_ID
    
    tmp.flj.from.snp.dat <- snp.gen.param.fun(snp.dat.indiv = tmp.snp.dat.indiv)
    
    rm(tmp.snp.dat.indiv)
    tmp.flj.from.snp.dat$FAM_SET_COMBN_ID <- fam.set.combn
    tmp.flj.from.snp.dat$SAMPLE_ID <- 0
    
    flj.from.snp.dat <- rbind(flj.from.snp.dat,tmp.flj.from.snp.dat)
  }
  
  rm(tmp.fams, tmp.flj.from.snp.dat, fam.set.combn)
  
  flj.from.snp.dat$A_TRANS_PROB <- flj.from.snp.dat$A_ALLELE_FREQ
  flj.from.snp.dat$B_TRANS_PROB <- flj.from.snp.dat$B_ALLELE_FREQ
  flj.from.snp.dat <- flj.from.snp.dat[,c("FAM_SET_COMBN_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB", "SAMPLE_ID")]
  
  return(flj.from.snp.dat)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

parent.combns.fun <- function(fams, 
                              n.in.pools,
                              fam.set.combns) {
  
  #generate parent.combns and fam.combns.by.fam.set.combn
  if(!is.null(fam.set.combns)) {
    fam.set.agg.combn <- unique(fam.set.combns[,c("FAM_SET_COMBN_ID", "FAM_SET_ID")])
    fam.set.agg.combn <- fam.set.agg.combn[order(fam.set.agg.combn[,"FAM_SET_ID"]),]
    fam.combns.by.fam.set.combn <- as.data.frame(matrix(NA, ncol = 1+n.in.pools, nrow=0))
    colnames(fam.combns.by.fam.set.combn) <- c("FAM_SET_COMBN_ID", paste("FAMILY_ID_", 1:n.in.pools, sep=""))
    
    for(fam.set.combn in unique(fam.set.combns[,"FAM_SET_COMBN_ID"])) {
      
      #pooling for pedigree assignment
      #  if(nrow(fam.set.agg.combn) == 
      #     n.in.pools * length(unique(fam.set.agg.combn$FAM_SET_COMBN_ID))) { 
      if(nrow(fam.set.agg.combn) != 1) {  
        for(agg in fam.set.agg.combn[fam.set.agg.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn,"FAM_SET_ID"]) {
          if(agg == min(fam.set.agg.combn[fam.set.agg.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn,"FAM_SET_ID"])) {
            tmp <- fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn &
                                    fam.set.combns[,"FAM_SET_ID"] == agg, 
                                  c("FAM_SET_COMBN_ID", "FAMILY_ID")]
          } else {
            tmp <- merge(tmp, fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn &
                                               fam.set.combns[,"FAM_SET_ID"] == agg, 
                                             c("FAM_SET_COMBN_ID", "FAMILY_ID")],
                         by = "FAM_SET_COMBN_ID",
                         all = TRUE,
                         suffixes = agg:(agg+ncol(tmp)))
          } 
        }
      } 
      
      #pooling by phenotype
      if(nrow(fam.set.agg.combn) == 1) { 
        for(agg in unique(fam.set.agg.combn$FAM_SET_ID)){  # 1:n.in.pools) {
          if(agg == unique(fam.set.agg.combn$FAM_SET_ID)[1]) {
            tmp <- fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn &
                                    fam.set.combns[,"FAM_SET_ID"] == agg, 
                                  c("FAM_SET_COMBN_ID", "FAMILY_ID")]
          } else {
            tmp <- merge(tmp, fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn , 
                                             c("FAM_SET_COMBN_ID", "FAMILY_ID")],
                         by = "FAM_SET_COMBN_ID",
                         all = TRUE,
                         suffixes = agg:(agg+ncol(tmp)))
          } 
        }
      }
      
      colnames(tmp) <- c("FAM_SET_COMBN_ID", paste("FAMILY_ID_", 1:n.in.pools, sep=""))
      fam.combns.by.fam.set.combn <- rbind(fam.combns.by.fam.set.combn, tmp)  
      rm(tmp)
    }
    
    parent.combns <- unique(fam.combns.by.fam.set.combn[,-1])
    
    if(!is.vector(parent.combns)) {
      parent.combns[,"FAM_COMBN_ID"] <- 1:nrow(parent.combns) 
    } else {
      parent.combns <- data.frame(FAMILY_ID_1 = parent.combns,
                                  FAM_COMBN_ID = 1:length(parent.combns)  )
      parent.combns$FAMILY_ID_1 <- as.integer(parent.combns$FAMILY_ID_1)
    }
    fam.combns.by.fam.set.combn <- merge(fam.combns.by.fam.set.combn, parent.combns, by.x = 2:ncol(fam.combns.by.fam.set.combn), 
                                         by.y = 1:(ncol(parent.combns)-1), all.x=TRUE)
    parent.combns <- parent.combns[,c("FAM_COMBN_ID", paste("FAMILY_ID_", 1:n.in.pools, sep=""))]
    
    fam.combns.by.fam.set.combn <- fam.combns.by.fam.set.combn[,c("FAM_SET_COMBN_ID", "FAM_COMBN_ID", paste("FAMILY_ID_", 1:n.in.pools, sep=""))]
    fam.combns.by.fam.set.combn <- fam.combns.by.fam.set.combn[order(fam.combns.by.fam.set.combn[, "FAM_COMBN_ID"]),]
    fam.combns.by.fam.set.combn <- fam.combns.by.fam.set.combn[order(fam.combns.by.fam.set.combn[, "FAM_SET_COMBN_ID"]),]
  }
  
  #Generate parent combinations and append to parent.combns and fam.combns.by.fam.set.combn
  
  for(col in 2:ncol(parent.combns)) {
    parent.combns <- merge(parent.combns, fams, by.x = col, by.y = "FAMILY_ID", all.x = TRUE, suffixes = (c(col-2,col-1)))
  }
  parent.combns[,(col+1):ncol(parent.combns)] <- t(apply(parent.combns[,(col+1):ncol(parent.combns)], 1, sort))
  colnames(parent.combns) <- c(colnames(parent.combns)[1:col],paste("PARENT_", 1:(2*(col-1)), sep=""))
  parent.combns <- parent.combns[order(parent.combns[,"FAM_COMBN_ID"]),]
  parent.combns[,"PARENT_COMBN_ID"] <- NA
  parent.combns[!duplicated(parent.combns[,(col+1):ncol(parent.combns)]) ,"PARENT_COMBN_ID"] <- 
    1:sum(!duplicated(parent.combns[,(col+1):ncol(parent.combns)]))
  #get unique parent combinations then remerge
  tmp <- parent.combns[!is.na(parent.combns[,"PARENT_COMBN_ID"]),(col+1):ncol(parent.combns)]
  parent.combns <- parent.combns[,colnames(parent.combns) != "PARENT_COMBN_ID"]
  parent.combns <- merge(parent.combns,
                         tmp,
                         by.x = (col+1):ncol(parent.combns),
                         by.y = 1:(ncol(tmp)-1),
                         all.x = TRUE)
  
  parent.combns <- parent.combns[,c("FAM_COMBN_ID",
                                    paste("FAMILY_ID_", 1:n.in.pools, sep =""),
                                    "PARENT_COMBN_ID",
                                    paste("PARENT_", 1:(2*n.in.pools), sep =""))]
  
  # parent.combns.by.fam.set.combn <- merge(fam.combns.by.fam.set.combn[,c("FAM_SET_COMBN_ID", "FAM_COMBN_ID")], 
  #                                 parent.combns[,c("FAM_COMBN_ID", "PARENT_COMBN_ID", paste("PARENT_", 1:(2*n.in.pools), sep =""))], 
  #                                 by = "FAM_COMBN_ID",
  #                                 all.x = TRUE)
  
  fam.combns.by.fam.set.combn$FAM_COMBN_ID <- as.integer(fam.combns.by.fam.set.combn$FAM_COMBN_ID)
  parent.combns$FAM_COMBN_ID      <- as.integer(parent.combns$FAM_COMBN_ID)
  parent.combns.by.fam.set.combn           <- left_join(fam.combns.by.fam.set.combn[,c("FAM_SET_COMBN_ID", "FAM_COMBN_ID")], 
                                                        parent.combns[,c("FAM_COMBN_ID", "PARENT_COMBN_ID", 
                                                                         paste("PARENT_", 1:(2*n.in.pools), sep =""))],
                                                        by = "FAM_COMBN_ID")
  
  parent.combns.by.fam.set.combn <- parent.combns.by.fam.set.combn[,colnames(parent.combns.by.fam.set.combn) != "FAM_COMBN_ID"]
  
  parent.combns.by.fam.set.combn <- parent.combns.by.fam.set.combn[!duplicated(parent.combns.by.fam.set.combn[,c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID")]),]
  
 parent.combns.by.fam.set.combn <- parent.combns.by.fam.set.combn[!duplicated(parent.combns.by.fam.set.combn[,c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID")]),] #remove duplicate parent combinations within family set combinations
 parent.combns <- parent.combns[order(parent.combns$FAM_COMBN_ID),]
 parent.combns <- parent.combns[!duplicated(parent.combns[,"PARENT_COMBN_ID"]),] #remove duplicate parent combinations 
  
  return(list(parent.combns = parent.combns, parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn))
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

lambda.ij.fun <- function(snp.dat.indiv, 
                          snp.param.indiv,
                          n.in.pools=1,
                          min.sd = 0,
                          min.intensity = 0
) {  
  
  #Returns the elements (e.g. "LAMBDA.AAAA", "LAMBDA.AAAB", "LAMBDA.AABB", "LAMBDA.BBBB") of the matrix shown in the 
  #top left of page 4 of Henshall et al. 2014 (modified to accomodate pools).  Note that it is assumed that 
  #LAMBDA.AAAB = LAMBDA.BAAA etc. Also returns the allelic proportion (pij) (see Equation 1 of Henshall et al. 2014). 
  
  #Required functions:
  # pij.fun
  
  # snp.dat.indiv: Data frame (pooled samples only) containing relevant fields from the Genotype.Intensity 
  #tab of the corresponding GenotypeIntensity.xls file outputted from Sequenom's Typer software 
  #(Sequenom 2006).  Equivalent outputs from other platforms could also be used.  
  #              1. SAMPLE_ID is the sample identifier, 
  #              2. SNP_ID   is the SNP identifier,
  #              3. INTENSITY_A   is the area/intensity for allele A, 
  #              4. INTENSITY_B   is the area/intensity for allele B, 
  
  # snp.param.indiv: Data frame. Output of snp.gen.param.fun.  Example below for n.in.pools = 2.
  #             SNP_ID	     SNP identifier
  #             MEAN_P_AAAA	 Mean of allelic proportion for genotype AAAA
  #             SD_P_AAAA	   Standard deviation of allelic proportion for genotype AAAA
  #             MEAN_P_AAAB	 Mean of allelic proportion for genotype AAAB
  #             SD_P_AAAB    Standard deviation of allelic proportion for genotype AAAB
  #             MEAN_P_AABB	 Mean of allelic proportion for genotype AABB
  #             SD_P_AABB	   Standard deviation of allelic proportion for genotype AABB
  #             MEAN_P_ABBB	 Mean of allelic proportion for genotype ABBB
  #             SD_P_ABBB    Standard deviation of allelic proportion for genotype ABBB
  #             MEAN_P_BBBB	 Mean of allelic proportion for genotype BBBB
  #             SD_P_BBBB	   Standard deviation of allelic proportion for genotype BBBB
  
  #n.in.pools
  # Integer. Number of individuals in each pooled sample.  NOTE computation time is roughly proportional to 2 n.in.pools (i.e. this function is only appropriate for small pool sizes).  For larger pool sizes use an alternative method (e.g. least squares).
  
  # min.sd       Number. Standard deviation of allelic proportion  fixed to this value if less 
  #              than it.
  
  # min.intensity      Number used in pij.fun. If sqrt((snp.dat.indiv$INTENSITY_A)^2 +
  #              (snp.dat.indiv$INTENSITY_B)^2) less than this value
  #              then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  #              Essentially removes observations that fall into the lower left of INTENSITY_A
  #              by INTENSITY_B scatter plot.
  
  #Returns##########################################
  # lambda.pool.j:  Data frame 
  #            1.            SAMPLE_ID is the individual identifier
  #            2.            SNP_ID is the SNP identifier 
  #            3-            e.g. AAAA_LAMBDA, AAAB_LAMBDA, AABB_LAMBDA, ABBB_LAMBDA, BBBB_LAMBDA - modified from top left 
  #                          of page 4 of Henshall et al. 2014 
  #            Final column. ALLELIC_PROP_POOL is the allellic proportion (Equation 1 of Henshall et al. 2014) 
  
  print("Running lambda.ij.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B") %in% 
         colnames(snp.dat.indiv)) != 4) {
    stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B")
  }
  
  #Generate list of genotypes
  genotypes <- genotypes.fun(n=(2*n.in.pools))
  
  if(sum(c("SNP_ID", paste("MEAN_P_",genotypes, sep=""), paste("SD_P_",genotypes, sep="")) %in% 
         colnames(snp.param.indiv)) != 3+4*n.in.pools) {
    stop(paste("snp.param.indiv input must be a data frame containing", 3+4*n.in.pools, "columns given that n.in.pools equals", n.in.pools))
  }
  
  #Name columns and assign class
  snp.dat.indiv$SAMPLE_ID  <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID   <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$INTENSITY_A   <- as.numeric(snp.dat.indiv$INTENSITY_A)
  snp.dat.indiv$INTENSITY_B   <- as.numeric(snp.dat.indiv$INTENSITY_B)
  snp.dat.indiv          <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")]
  
  #colnames(snp.param.indiv) <- c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB")
  #  snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  #  snp.param.indiv$MEAN_P_AA   <- as.numeric(snp.param.indiv$MEAN_P_AA)
  #  snp.param.indiv$SD_P_AA     <- as.numeric(snp.param.indiv$SD_P_AA)
  #  snp.param.indiv$MEAN_P_AB   <- as.numeric(snp.param.indiv$MEAN_P_AB)
  #  snp.param.indiv$SD_P_AB     <- as.numeric(snp.param.indiv$SD_P_AB)
  #  snp.param.indiv$MEAN_P_BB   <- as.numeric(snp.param.indiv$MEAN_P_BB)
  #  snp.param.indiv$SD_P_BB     <- as.numeric(snp.param.indiv$SD_P_BB)
  #  snp.param.indiv <- snp.param.indiv[,c("SNP_ID", "MEAN_P_AA", "SD_P_AA", "MEAN_P_AB", "SD_P_AB", "MEAN_P_BB", "SD_P_BB")]
  
  #Check for duplicated records in snp.dat.indiv
  pool.snp <- paste(snp.dat.indiv$SAMPLE_ID,snp.dat.indiv$SNP_ID, sep=".")
  if(sum(duplicated(pool.snp)) > 0) {
    stop("SAMPLE_ID and SNP_ID combinations are not unique in snp.dat.indiv.  Delete duplicates or recode SAMPLE_ID.")
  }
  rm(pool.snp)
  
  # Check the list of SNPs the same in input files
  if(sum(unique(snp.dat.indiv[,"SNP_ID"]) != unique(snp.param.indiv[,"SNP_ID"]))>0) {
    stop("SNP identifiers do not match in snp.dat.indiv and snp.param.indiv")
  }
  
  #Get allelic proportion
  snp.allelic.prop <- pij.fun(snp.dat.indiv = snp.dat.indiv, min.intensity = min.intensity)
  print("Still running lambda.ij.fun")
  
  #Rename columns
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "SAMPLE_ID"]        <- "SAMPLE_ID"
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "ALLELIC_PROP"]     <- "ALLELIC_PROP_POOL"  
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "INTENSITY"]        <- "INTENSITY_POOL"  
  
  #merge data using cbind
  if(identical(snp.dat.indiv$SAMPLE_ID,snp.allelic.prop$SAMPLE_ID) &
     identical(snp.dat.indiv$SNP_ID,snp.allelic.prop$SNP_ID)) {
    snp.dat.indiv <- cbind(snp.dat.indiv, snp.allelic.prop[,!colnames(snp.allelic.prop) %in% c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")])
  } else {
    stop("SAMPLE_ID and SNP_ID columns of pij.fun output do not match those of snp.dat.indiv.  Not sure why.")
  }
  
  for(colname in paste("SD_P_",genotypes, sep="")) {
    snp.param.indiv[snp.param.indiv[,colname] < min.sd & !is.na(snp.param.indiv[,colname]), colname] <- min.sd
    snp.param.indiv[is.na(snp.param.indiv[,colname]), colname] <- min.sd
  }
  
  for(genotype in genotypes) {
    mean.colname <- paste("MEAN_P_",genotype, sep="")
    sd.colname <- paste("SD_P_",genotype, sep="")
    
    #Ensure all SD are > 0 if mean is not NA
    if(sum(
      (
        !is.na(snp.param.indiv[,mean.colname]) &
        (snp.param.indiv[,sd.colname] <= 0 & !is.na(snp.param.indiv[,sd.colname]) |
         is.na(snp.param.indiv[,sd.colname])) 
      )
    ) > 0) {
      stop("At least one element of SD_P_ equals zero or is NA where MEAN_P_ are not NA in snp.param.indiv.  Must specify min.sd that is greater than 0 or modify snp.param.indiv.")
    }
  }
  rm(mean.colname, sd.colname)
  
  #Merge data
  # lambda.pools <- merge(snp.dat.indiv, snp.param.indiv, by = "SNP_ID", all.x = TRUE)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  lambda.pools <- left_join(snp.dat.indiv, snp.param.indiv, by = "SNP_ID")
  
  for(count in 1:length(genotypes)) {
    genotype     <- genotypes[count]
    lambda.colname  <- paste(genotype, "_LAMBDA",sep="")
    mean.colname <- paste("MEAN_P_",genotype, sep="")
    sd.colname   <- paste("SD_P_",genotype, sep="")
    
    lambda.pools[,lambda.colname] <- dnorm(lambda.pools$ALLELIC_PROP_POOL,lambda.pools[,mean.colname],lambda.pools[,sd.colname])
  }
  
  #retain only relevant columns
  lambda.pools <- lambda.pools[,c("SAMPLE_ID", "SNP_ID", paste(genotypes, "_LAMBDA",sep=""), "ALLELIC_PROP_POOL")]
  
  return(lambda.pools)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

rho.inv.fun <- function(n.in.pools) {
  
  print("Running rho.inv.fun")
  
  #Generate list of genotypes
  genotypes <- genotypes.fun(n=(2*n.in.pools))
  
  rho.inv <- data.frame(GENOTYPE = genotypes,
                        RHO_INV      = NA)
  rho.inv$GENOTYPE <- as.character(rho.inv$GENOTYPE)
  
  combns.fun <- function(n, m) {
    ind <- combn(seq_len(n), m)
    ind <- t(ind) + (seq_len(ncol(ind)) - 1) * n
    res <- rep(0, nrow(ind) * n)
    res[ind] <- 1
    matrix(res, ncol = n, nrow = nrow(ind), byrow = TRUE)
  }
  
  for(i in 1:length(genotypes)) {
    geno <- genotypes[i]
    ordered.geno.count <- nrow(combns.fun(2*n.in.pools,lengths(regmatches(geno, gregexpr("B", geno)))))
    rho.inv[i,"RHO_INV"] <- 1/ordered.geno.count
    rm(ordered.geno.count)
  }
  
  return(rho.inv)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

nlj.fun <- function(flj, n.in.pools) {
  print("Running nlj.fun")
  
  #Get nlj
  tmp.parent.combns <- data.frame(PARENT_COMBN_ID = 1,
                                  PARENT_1 = 1)
  
  for(i in 2:(2*n.in.pools)) {
    tmp.parent.combns[,paste("PARENT_", i, sep="")] <- 1
  }  
  
  nlj <- NULL
  for(fam.set.combn in unique(flj$FAM_SET_COMBN_ID)) {
    
    tmp.fj <- flj[flj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    
    tmp.tij <- tmp.fj
    tmp.tij$SAMPLE_ID = 1
    
    tmp.nj <- tcj.fun(parent.combns = tmp.parent.combns,
                      tij = tmp.tij,
                      fj = tmp.fj) 
    tmp.nj <- tmp.nj[,!colnames(tmp.nj) %in% c("PARENT_COMBN_ID", "MISS_PARENT_COUNT")]
    tmp.nj$FAM_SET_COMBN_ID <- fam.set.combn
    
    nlj <- rbind(nlj,tmp.nj)
  }
  
  rm(tmp.fj, tmp.parent.combns, tmp.tij, tmp.nj)
  
  #reorder columns
  nlj <- nlj[,c(colnames(nlj)[ncol(nlj)], colnames(nlj)[-ncol(nlj)])]
  
  print("Still running nlj.fun")
  
  # }
  
  return(nlj)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

tcj.fun <- function(tij, parent.combns, fj = NULL) {
  
  #Args##########################################
  
  #tij: Data frame.  Output of Gij.fun,  Mij.fun or Dij.fun.
  #          1.   SAMPLE_ID is the individual identifier (assumed to be offspring if not 
  #               listed as a sire or dam in 'fams' input)
  #          2.   SNP_ID is the SNP identifier
  #          3-4. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #               for alleles A and B respectively computed from Gij (or Mij) (i.e. the elements of 
  #               the transmission (Tij) vector, Equation 2 of Henshall et al. 2014).
  #          5.    OFFSPRING_MISS. Logical. True if genotype data are missing
  
  #NOTE: If a D.matrix is specified then fj must be NULL (fj = NULL)
  
  # parent.combns:   
  
  # fj: Data frame.  Output of fj.from.parent.Gij.fun or fj.from.snp.dat.indiv.fun
  #          1.   SNP_ID is the snp identifier
  #          2-3. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #               for alleles A and B respectively computed from average parental Gij or parental genotype counts (i.e. the elements of 
  #               the fj vector).
  #          4. SAMPLE_ID is the individual identifier (all values are zero and represent unknown individuals)
  
  #Returns##########################################
  
  #tcj: Data frame.
  
  
  print("Running tcj.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB") %in% 
         colnames(tij)) != 4) {
    stop("tij input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, A_TRANS_PROB, B_TRANS_PROB")
  }
  
  #Name columns and assign class
  #colnames(tij)     <- c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")
  tij$SAMPLE_ID     <- as.integer(tij$SAMPLE_ID)
  tij$SNP_ID       <- as.character(tij$SNP_ID)
  tij$A_TRANS_PROB <- as.numeric(tij$A_TRANS_PROB)
  tij$B_TRANS_PROB <- as.numeric(tij$B_TRANS_PROB)
  tij <- tij[,c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")]
  
  if(is.null(fj)) {
    fj <- data.frame(SNP_ID       = unique(tij$SNP_ID),
                     A_TRANS_PROB = 0.5,
                     B_TRANS_PROB = 0.5) #Missing parent data for exclusion method
  } else {
    fj$SAMPLE_ID     <- as.integer(fj$SAMPLE_ID)
    fj$SNP_ID       <- as.character(fj$SNP_ID)               
    fj$A_TRANS_PROB <- as.numeric(fj$A_TRANS_PROB)  
    fj$B_TRANS_PROB <- as.numeric(fj$B_TRANS_PROB)  
    fj <- fj[,c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")]
    
  }
  
  #Name columns and assign class
  for(i in 1:ncol(parent.combns)) {parent.combns[,i] <- as.integer(parent.combns[,i])}
  parent.combns$PARENT_COMBN_ID <- as.integer(parent.combns$PARENT_COMBN_ID)
  
  #Check that probabilities between 0 and 1
  if(
    sum(
      tij$A_TRANS_PROB > 1 |
      tij$B_TRANS_PROB > 1 |
      
      tij$A_TRANS_PROB < 0 |
      tij$B_TRANS_PROB < 0 
      , na.rm = T) != 0
  ) {
    stop("Probabilities in tij must be between 0 and 1 inclusive")
  }
  
  #Check that "A_TRANS_PROB", "B_TRANS_PROB" sum to 1 in tij
  if(
    sum(
      round((tij$A_TRANS_PROB + tij$B_TRANS_PROB),5) != 1.0, na.rm = T
    ) != 0
  ) {
    stop("A_TRANS_PROB + B_TRANS_PROB must equal 1 in all rows of tij")
  }
  
  #Ensure that SAMPLE_ID column contains no 0s or NA
  if(
    sum(
      (is.na(tij$SAMPLE_ID) |
       tij$SAMPLE_ID == 0), na.rm = T
    ) != 0
  ) {
    stop("SAMPLE_ID in tij cannot be NA or 0")
  } 
  
  #remove family combinations from parent.combns
  parent.combns <- parent.combns[!duplicated(parent.combns[,"PARENT_COMBN_ID"]),]
  parent.combns <- parent.combns[,grepl("PARENT_", colnames(parent.combns))]
  
  #identify parents 
  parents <- unique(unlist(parent.combns[,colnames(parent.combns) != "PARENT_COMBN_ID"]))
  
  #separate parents and offspring in tij
  tij.parents <- tij[tij[,"SAMPLE_ID"] %in% parents,]
  
  #Include SAMPLE_ID = 0 (i.e. unknown parent)
  if(!is.null(fj)) {
    tij.parents      <- merge(tij.parents, fj, 
                              by = c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB"), 
                              all = TRUE)
  }
  
  #identify snp
  snps <- fj$SNP_ID
  
  #Transmission vector 
  tcj <- parent.combns[rep(seq_len(nrow(parent.combns)), length(snps)), ] 
  tcj[,"SNP_ID"] <- rep(snps, each = nrow(parent.combns))
  rownames(tcj) <- 1:nrow(tcj)
  gc()
  
  tcj$MISS_PARENT_COUNT <- 0
  tcj$MISS <- FALSE
  
  # tcj <- merge(tcj, fj[,c("SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")], by = "SNP_ID", all.x = TRUE)
  tcj <- left_join(tcj, fj[,c("SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")], by = "SNP_ID")
  
  colnames(tcj)[(ncol(tcj)-1):ncol(tcj)] <- c("A_TRANS_PROB_MISS", "B_TRANS_PROB_MISS")
  
  for(parent in 1:(ncol(parent.combns) - 1)) {
    
    print(paste("Parent",parent,"of",(ncol(parent.combns) - 1)))
    
    # tcj <- merge(tcj,
    #             tij.parents[,c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")],
    #             by.x = c(paste("PARENT_", parent, sep=""), "SNP_ID"),
    #              by.y = c("SAMPLE_ID", "SNP_ID"),
    #              all.x = TRUE)
    
    tmp <- c("PARENT" = "SAMPLE_ID", "SNP_ID" = "SNP_ID")
    attributes(tmp)$names[1] <- paste("PARENT_", parent, sep="")
    tcj <- left_join(tcj, tij.parents[,c("SAMPLE_ID", "SNP_ID", "A_TRANS_PROB", "B_TRANS_PROB")], by = tmp)
    rm(tmp)
    
    tcj$MISS <- FALSE
    tcj[is.na(tcj[,"A_TRANS_PROB"]) | 
          is.na(tcj[,"B_TRANS_PROB"]) |
          tcj[,paste("PARENT_", parent, sep="")] == 0,"MISS"] <- TRUE
    
    tcj$MISS_PARENT_COUNT <- tcj$MISS_PARENT_COUNT + tcj$MISS
    
    #Replace missing data with mean allele probabilities in parents
    tcj[is.na(tcj[,"A_TRANS_PROB"]),"A_TRANS_PROB"] <-
      tcj[is.na(tcj[,"A_TRANS_PROB"]),"A_TRANS_PROB_MISS"]
    tcj[is.na(tcj[,"B_TRANS_PROB"]),"B_TRANS_PROB"] <-
      tcj[is.na(tcj[,"B_TRANS_PROB"]),"B_TRANS_PROB_MISS"]
    
    if(parent == 1) {
      colnames(tcj)[(ncol(tcj)-1):ncol(tcj)] <- c("A","B") #last 2 columns
    } else {
      tcj[,paste(genotypes.fun(n = (parent-1)),"A", sep="")] <- 
        tcj[,genotypes.fun(n = (parent-1))] * tcj[,"A_TRANS_PROB"]
      tcj[,paste(genotypes.fun(n = (parent-1)),"B", sep="")] <- 
        tcj[,genotypes.fun(n = (parent-1))] * tcj[,"B_TRANS_PROB"]
      
      #remove unwanted columns
      tcj <- tcj[,!colnames(tcj) %in% c("A_TRANS_PROB", "B_TRANS_PROB")]
      tcj <- tcj[,!colnames(tcj) %in% genotypes.fun(parent-1)]
      
      #convert ordered genotypes to unordered
      ordered.genos <- c(paste(genotypes.fun(n = (parent-1)),"A", sep=""),
                         paste(genotypes.fun(n = (parent-1)),"B", sep=""))
      unordered.genos <- genotypes.fun(n = parent)
      
      for (i in 1:length(ordered.genos)) {
        unordered.genos[i] <- paste(sort(unlist(strsplit(ordered.genos[i], ""))), collapse = "")
      }
      
      #partition tcj 
      tcj.ordered <- tcj[,ordered.genos]
      tcj <- tcj[,!colnames(tcj) %in% ordered.genos]
      
      #get retain unordered genotypes
      for(geno in unique(unordered.genos)) {
        if(sum(unordered.genos==geno) > 1) {
          tcj[,geno] <- rowSums(tcj.ordered[,unordered.genos==geno])
        } else {
          tcj[,geno] <- tcj.ordered[,unordered.genos==geno]
        }
      }
      rm(tcj.ordered)
    }
    
  }
  #reorder/remove columns
  tcj <- tcj[,c("SNP_ID", "PARENT_COMBN_ID", "MISS_PARENT_COUNT", genotypes.fun(parent))]
  
  #divide by rho as genotypes are unordered
  n.in.pools <- (ncol(parent.combns)-1)/2
  
  rho.inv <- rho.inv.fun(n.in.pools)
  for(geno in rho.inv$GENOTYPE) {
    tcj[,geno] <- tcj[,geno]*rho.inv[rho.inv[,"GENOTYPE"] == geno,"RHO_INV"]
  }
  
  tcj <- tcj[order(tcj[,"SNP_ID"]),]
  tcj <- tcj[order(tcj[,"PARENT_COMBN_ID"]),]
  
  return(tcj)
  
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

snp.error.fun <- function(Gij, fams = NULL, parents.only = FALSE) {    
  
  #Computes the estimated SNP specific error rate (epsilon hat) from Gij matrices according to the method
  #outlined on the left hand side of  page 5 of Henshall et al. 2014 
  
  # Args: 
  
  # Gij:     Data frame.  Output from G.T.matrices.fun
  #              1.   SAMPLE_ID is the individual identifier
  #              2.   SNP_ID is the SNP identifier
  #              3-6. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the Gij matrix (see the 
  #                   top left of page 4 of Henshall et al. 2014
  #              7-8. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #                   for alleles A and B respectively computed from Gij (i.e. the elements of 
  #                   the transmission (Tij) vector, Equation 2 of Henshall et al. 2014).
  
  # fams: Optional data frame - only required if parents.only=TRUE   NOTE: to define half sib families, the unknown parent must be NA or 0.  
  #          1. FAMILY_ID is the family identifier
  #          2. SIRE_ID is the sire identifier 
  #          3. DAM_ID is the dam identifier 
  
  # parents.only: Logical.  If TRUE only used data from parents.
  
  #Returns##########################################
  
  # snp.error.hat: Data frame.
  #              1.   SNP_ID is the SNP identifier  
  #              2.   SNP_ERROR_HAT is the SNP specific error rate (epsilon hat) from 
  #                   Gij matrix according to the method outlined on the left hand side of  
  #                   page 5 of Henshall et al. 2014  
  
  print("Running snp.error.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", 
           "BB_GENO_PROB", "A_TRANS_PROB", "B_TRANS_PROB") %in% 
         colnames(Gij)) != 8) {
    stop("Gij input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB, BB_GENO_PROB, A_TRANS_PROB, B_TRANS_PROB")
  }
  
  
  #Name columns and assign class
  #colnames(Gij)     <- c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB", 
  #                                 "A_TRANS_PROB", "B_TRANS_PROB")
  Gij$SAMPLE_ID     <- as.integer(Gij$SAMPLE_ID)
  Gij$SNP_ID       <- as.character(Gij$SNP_ID)
  Gij$AA_GENO_PROB <- as.numeric(Gij$AA_GENO_PROB)
  Gij$AB_GENO_PROB <- as.numeric(Gij$AB_GENO_PROB)
  Gij$BA_GENO_PROB <- as.numeric(Gij$BA_GENO_PROB)
  Gij$BB_GENO_PROB <- as.numeric(Gij$BB_GENO_PROB)
  Gij$A_TRANS_PROB <- as.numeric(Gij$A_TRANS_PROB)
  Gij$B_TRANS_PROB <- as.numeric(Gij$B_TRANS_PROB)
  Gij <- Gij[,c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", 
                "BB_GENO_PROB", "A_TRANS_PROB", "B_TRANS_PROB")]
  
  #Check that probabilities between 0 and 1
  if(
    sum(
      Gij$AA_GENO_PROB > 1 |
      Gij$AB_GENO_PROB > 1 |
      Gij$BA_GENO_PROB > 1 |
      Gij$BB_GENO_PROB > 1 |
      Gij$A_TRANS_PROB > 1 |
      Gij$B_TRANS_PROB > 1 |
      
      Gij$AA_GENO_PROB < 0 |
      Gij$AB_GENO_PROB < 0 |
      Gij$BA_GENO_PROB < 0 |
      Gij$BB_GENO_PROB < 0 |
      Gij$A_TRANS_PROB < 0 |
      Gij$B_TRANS_PROB < 0 
      , na.rm = T) != 0
  ) {
    stop("Probabilities in Gij must be between 0 and 1 inclusive")
  }
  
  #Check that "A_TRANS_PROB", "B_TRANS_PROB" sum to 1 in Gij
  if(
    sum(
      round((Gij$A_TRANS_PROB + Gij$B_TRANS_PROB),5) != 1.0, na.rm = T
    ) != 0
  ) {
    stop("A_TRANS_PROB + B_TRANS_PROB must equal 1 in all rows of Gij")
  }
  
  #Check that "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB" sum to one
  if(
    sum(
      round((Gij$AA_GENO_PROB + Gij$AB_GENO_PROB +
             Gij$BA_GENO_PROB + Gij$BB_GENO_PROB),5) != 1, na.rm = T
    ) != 0
  ) {
    stop("AA_GENO_PROB + AB_GENO_PROB + BA_GENO_PROB + BB_GENO_PROB must equal 1 in all rows of Gij")
  }
  
  #Ensure that SAMPLE_ID column contains no 0s or NA
  if(
    sum(
      (is.na(Gij$SAMPLE_ID) |
       Gij$SAMPLE_ID == 0), na.rm = T
    ) != 0
  ) {
    stop("SAMPLE_ID in Gij cannot be NA or 0")
  } 
  
  #only retain parents in Gij if parents.only = TRUE
  
  if(parents.only == TRUE) {
    
    if(is.null(fams)) {
      stop("fams must be specified if parents.only = TRUE")
    }
    
    if(sum(c("FAMILY_ID", "SIRE_ID", "DAM_ID") %in% 
           colnames(fams)) != 3) {
      stop("fams input must be a data frame containing the following headings: FAMILY_ID, SIRE_ID, DAM_ID")
    }
    
    #colnames(fams) <- c("FAMILY_ID", "SIRE_ID", "DAM_ID")
    fams$FAMILY_ID    <- as.integer(fams$FAMILY_ID)
    fams$SIRE_ID   <- as.integer(fams$SIRE_ID)
    fams$DAM_ID    <- as.integer(fams$DAM_ID)
    fams <- fams[,c("FAMILY_ID", "SIRE_ID", "DAM_ID")]
    
    #If NA present as parents in fam then convert to 0
    fams[is.na(fams[,"SIRE_ID"]) ,"SIRE_ID"] <- "0"   
    fams[is.na(fams[,"DAM_ID"]) ,"DAM_ID"]   <- "0"  
    
    parents <- unique(c(fams[,"SIRE_ID"],fams[,"DAM_ID"]))
    parents <- parents[parents != "0"]  
    
    Gij <- Gij[Gij[,"SAMPLE_ID"] %in% parents,]
    
  }
  
  #Compute 1 - maximum of the unordered genotypes in Gij
  Gij$ab.BA_GENO_PROB   <- Gij$AB_GENO_PROB + Gij$BA_GENO_PROB
  Gij$EPSILON_GENO_PROB <- 1 - do.call(pmax, Gij[,c("AA_GENO_PROB", 
                                                    "ab.BA_GENO_PROB", 
                                                    "BB_GENO_PROB")])
  #Compute SNP error rates (mean of EPSILON_GENO_PROB for each SNP)
  snp.error <- aggregate(Gij$EPSILON_GENO_PROB, 
                         by = list(Gij$SNP_ID), 
                         na.rm=T, 
                         FUN = "mean")      
  colnames(snp.error) <- c("SNP_ID", "SNP_ERROR_HAT")
  
  return(snp.error)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

prelim.ml.quant.fun <- function(Gij,
                                parent.combns.by.fam.set.combn,
                                flj,
                                gkj,
                                snp.error,
                                nlj,
                                fam.set.combns.by.pool) {
  
  print("Running prelim.ml.quant.fun")
  
  tclj <- NULL
  for(fam.set.combn in unique(flj$FAM_SET_COMBN_ID)) {
    
    print(fam.set.combn)
    
    tmp.fj <- flj[flj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    tmp.parent.combns <- parent.combns.by.fam.set.combn[parent.combns.by.fam.set.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn,] 
    
    tmp.parents <- tmp.parent.combns[,!colnames(tmp.parent.combns) %in% c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID")]
    tmp.parents <- unique(as.vector(as.matrix(tmp.parents)))
    
    tmp.tij <- Gij[Gij[,"SAMPLE_ID"] %in% tmp.parents, ]
    tmp.tclj <- tcj.fun(tij = tmp.tij, 
                        parent.combns = tmp.parent.combns, 
                        fj = tmp.fj)
    tmp.tclj$FAM_SET_COMBN_ID <- fam.set.combn
    tclj <- rbind(tclj, tmp.tclj)
  }
  rm(tmp.fj, tmp.parent.combns, tmp.parents, tmp.tij, tmp.tclj)
  
  miss.parent.count <- unique(tclj[,c("PARENT_COMBN_ID", "SNP_ID", "MISS_PARENT_COUNT")])
  tclj               <- tclj[,colnames(tclj) != "MISS_PARENT_COUNT"] 
  
  gklj   <- adj.geno.prob.fun(lhs = merge(gkj,fam.set.combns.by.pool, by = "SAMPLE_ID", all = TRUE),
                              snp.error = snp.error,
                              nlj = nlj)
  return(list(tclj = tclj,
              miss.parent.count = miss.parent.count,
              gklj = gklj))
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

adj.geno.prob.fun <- function(lhs, snp.error, nlj) {
  
  #Computes the Mij matrix (see the top left of page 5 of Henshall et al. 2014)
  
  #Args########################################## 
  #lhs. Data frame.  
  #              1.   FAM_SET_COMBN_ID
  #              2.   SAMPLE_ID 
  #              3.   SNP_ID is the SNP identifier
  #              4-.  Unsorted genotypes: elements of the gkj or tclj matrices 
  
  #snp.error. Data frame. Output from snp.error.fun
  #              1.   SNP_ID is the SNP identifier
  #              2.   SNP_ERROR_TILDE is the assumed SNP specific error rate (see the 
  #                   middle left of page 5 of Henshall et al. 2014)
  
  #nlj.  Data frame.  Output of nlj.fun
  #              1.   SNP_ID is the SNP identifier
  #              2-   Unsorted genetypes
  
  #Returns##########################################
  
  # adj.geno.probs:     Data frame 
  #              1.   SAMPLE_ID or FAM_COMBN_ID 
  #              2.   SNP_ID is the SNP identifier
  #              3-.  Unsorted genotypes: elements of the gkj or tclj matrices 
  
  print("Running adj.geno.prob.fun")
  
  #Name columns and assign class
  
  # if(sum(c("SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB") %in% colnames(fjfj)) != 5) {
  #    stop("fjfj input must be a data frame containing the following headings: SNP_ID, AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB, BB_GENO_PROB")
  #  }
  
  #  fjfj$SNP_ID       <- as.character(fjfj$SNP_ID)
  #  fjfj$AA_GENO_PROB <- as.numeric(fjfj$AA_GENO_PROB)
  #  fjfj$AB_GENO_PROB <- as.numeric(fjfj$AB_GENO_PROB)
  #  fjfj$BA_GENO_PROB <- as.numeric(fjfj$BA_GENO_PROB)
  #  fjfj$BB_GENO_PROB <- as.numeric(fjfj$BB_GENO_PROB)
  #  fjfj <- fjfj[,c("SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")]

  if(sum(c("SNP_ID", "SNP_ERROR_TILDE") %in% colnames(snp.error)) != 2) {
    stop("snp.error input must be a data frame containing the following headings: SNP_ID, SNP_ERROR_TILDE")
  }
  
  snp.error$SNP_ID       <- as.character(snp.error$SNP_ID)
  snp.error$SNP_ERROR_TILDE <- as.numeric(snp.error$SNP_ERROR_TILDE)
  snp.error <- snp.error[,c("SNP_ID", "SNP_ERROR_TILDE")]
  
  #  if(sum(c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB") %in% colnames(lhs)) != 6) {
  #    stop("lhs input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB, BB_GENO_PROB")
  #  }
  
  #  lhs$SAMPLE_ID     <- as.integer(lhs$SAMPLE_ID)
  #  lhs$SNP_ID       <- as.character(lhs$SNP_ID)
  #  lhs$AA_GENO_PROB <- as.numeric(lhs$AA_GENO_PROB)
  #  lhs$AB_GENO_PROB <- as.numeric(lhs$AB_GENO_PROB)
  #  lhs$BA_GENO_PROB <- as.numeric(lhs$BA_GENO_PROB)
  #  lhs$BB_GENO_PROB <- as.numeric(lhs$BB_GENO_PROB)
  #  lhs <- lhs[,c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")]
  nlj <- nlj[order(nlj$SNP_ID),]
  lhs <- lhs[order(lhs$SNP_ID),]
  snp.error <- snp.error[order(snp.error$SNP_ID),]
  
  #Check that SNP_IDs the same in all data frames
  if(!(identical(unique(nlj$SNP_ID), unique(lhs$SNP_ID)))) {
    stop("SNP_ID in nlj must match those in lhs")
  }
  
  if(!(identical(unique(snp.error$SNP_ID), unique(lhs$SNP_ID)))) {
    stop("SNP_ID in snp.error must match those in lhs")
  }
  
  #Check that SNP_IDs the same in all data frames
  if(!(identical(unique(nlj$SNP_ID), unique(snp.error$SNP_ID)))) {
    stop("SNP_ID in nlj must match those in snp.error")
  }
  
  #End checks  
  
  # remove from consideration snp that have a 0 in any column of nj 
  snp.remove <- nlj[rowSums(nlj[,!colnames(nlj) %in% c("FAM_SET_COMBN_ID", "SNP_ID")] == 0) > 0,c("FAM_SET_COMBN_ID", "SNP_ID")]
  
  if(length(snp.remove) > 0) {
    for(set in unique(snp.remove[,"FAM_SET_COMBN_ID"])) {
      
      nlj <- nlj[!nlj[,"SNP_ID"] %in% snp.remove[snp.remove[,"FAM_SET_COMBN_ID"] == set,"SNP_ID"] & 
                   nlj[,"FAM_SET_COMBN_ID"] == set,]
      lhs <- lhs[!lhs[,"SNP_ID"] %in% snp.remove[snp.remove[,"FAM_SET_COMBN_ID"] == set,"SNP_ID"] & 
                   lhs[,"FAM_SET_COMBN_ID"] == set,]
      warning(paste("Removed SNP",snp.remove[snp.remove[,"FAM_SET_COMBN_ID"] == set,"SNP_ID"], "from tcj.adj and gkj.adj of FAM_SET_COMBN_ID", set, "as zeros present in nlj)"))
    }
  }
  
  adj.geno.probs <- NULL
  for (fam.set.combn in unique(nlj$FAM_SET_COMBN_ID)) {
    
    tmp.lhs <- lhs[lhs[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    tmp.lhs <- tmp.lhs[,colnames(tmp.lhs) != "FAM_SET_COMBN_ID"]
    tmp.nlj <- nlj[nlj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]   
    tmp.nlj <- tmp.nlj[,colnames(tmp.nlj) != "FAM_SET_COMBN_ID"]  
    
    #Replace missing values in tmp.lhs with those from tmp.nlj
    tmp.lhs <- merge(tmp.lhs, tmp.nlj, by = "SNP_ID", suffix = c("","y"), all.x = TRUE) 
    tmp.lhs[,3:(ncol(tmp.lhs)/2+1)][rowSums(is.na(tmp.lhs[,3:(ncol(tmp.lhs)/2+1)])) > 0,] <- 
      tmp.lhs[,(ncol(tmp.lhs)/2+2):ncol(tmp.lhs)][rowSums(is.na(tmp.lhs[,3:(ncol(tmp.lhs)/2+1)])) > 0,]
    
    #Get SNP error
    # tmp.lhs <- merge(tmp.lhs, snp.error, by = "SNP_ID", all.x = TRUE)
    tmp.lhs$SNP_ID       <- as.character(tmp.lhs$SNP_ID)
    snp.error$SNP_ID <- as.character(snp.error$SNP_ID)
    tmp.lhs <- left_join(tmp.lhs, snp.error, by = "SNP_ID")
    
    #Make adjustments to gkj and nj to account for SNP error 
    tmp.lhs[,3:((ncol(tmp.lhs)-1)/2+1)] <- tmp.lhs[,3:((ncol(tmp.lhs)-1)/2+1)]*(1-tmp.lhs[,ncol(tmp.lhs)])
    tmp.lhs[,((ncol(tmp.lhs)-1)/2+2):(ncol(tmp.lhs)-1)]   <- tmp.lhs[,((ncol(tmp.lhs)-1)/2+2):(ncol(tmp.lhs)-1)]*tmp.lhs[,ncol(tmp.lhs)] 
    tmp.lhs <- tmp.lhs[,-(ncol(tmp.lhs))] #remove error column
    tmp.lhs[,3:(ncol(tmp.lhs)/2+1)] <- tmp.lhs[,3:(ncol(tmp.lhs)/2+1)] + tmp.lhs[,(ncol(tmp.lhs)/2+2):ncol(tmp.lhs)]
    tmp.adj.geno.probs <- tmp.lhs[,1:(ncol(tmp.lhs)/2+1)]
    
    tmp.adj.geno.probs <- tmp.adj.geno.probs[order(tmp.adj.geno.probs[,"SNP_ID"]),]
    
    tmp.adj.geno.probs$FAM_SET_COMBN_ID <- fam.set.combn
    adj.geno.probs <- rbind(adj.geno.probs,tmp.adj.geno.probs)
  }
  
  return(adj.geno.probs)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

ml.fun <- function(g.d.klj.adj,
                   tclj,
                   snp.error,
                   nlj,
                   miss.parent.count,
                   parent.combns.by.fam.set.combn,
                   parent.combns,
                   meth) {
  
  print("Running ml.fun")
  
  wd <- getwd()
  
  dir.create(file.path(wd, "Results"), showWarnings = FALSE)
  setwd(file.path(wd, "Results"))
  
  g.d.klj.adj <- g.d.klj.adj[order(g.d.klj.adj[,"SAMPLE_ID"]),]
  
  tclj.adj <- adj.geno.prob.fun(lhs= tclj,
                                snp.error=snp.error,
                                nlj = nlj)
  
  tclj.adj <- tclj.adj[order(tclj.adj[,"PARENT_COMBN_ID"]),]
  
  duos <- parents.to.pools.lod.fun(g.d.klj.adj = g.d.klj.adj,
                                   tclj.adj = tclj.adj,  
                                   miss.parent.count = miss.parent.count,
                                   nlj = nlj,
                                   parent.combns.by.fam.set.combn = parent.combns.by.fam.set.combn)
  duos.lod  <- duos$duos.lod
  duos.logl <- duos$duos.logl
  rm(duos)
  
  # duos.lod <- merge(duos.lod, parent.combns.by.fam.set.combn, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.x = TRUE)
  duos.lod <- left_join(duos.lod, unique(g.d.klj.adj[,c("SAMPLE_ID", "FAM_SET_COMBN_ID")]), by = "SAMPLE_ID")
  duos.lod <- left_join(duos.lod, parent.combns.by.fam.set.combn, by = c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID"))
  
  duos.lod <- duos.lod[order(duos.lod$LOD, decreasing = TRUE),]
  duos.lod <- duos.lod[order(duos.lod$SAMPLE_ID),]
  
  most.like.parents <- most.like.parents.duo.fun(duos.lod = duos.lod,
                                                 parent.combns = parent.combns) 
  
  if(!exists("running.sim")) {
    try(plot.lod.vs.LogL.range.fun(lod = duos.lod[,c("SAMPLE_ID", "RANGE_5_TO_95_LOGL", "LOD")],
                                   meth = meth)
    )
  }
  setwd(wd)
  
  return(list(g.d.klj.adj = g.d.klj.adj,
              tclj.adj = tclj.adj,
              duos.lod = duos.lod,
              duos.logl = duos.logl,
              most.like.parents = most.like.parents))
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

parents.to.pools.lod.fun <- function(g.d.klj.adj,
                                     tclj.adj,  
                                     miss.parent.count = NULL,
                                     nlj,
                                     parent.combns.by.fam.set.combn) {  
  
  # fjfj: Data frame.  Output of flj.from.parent.Gij.fun or flj.from.snp.dat.fun
  #          1.   SNP_ID is the snp identifier
  #          2-3. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #               for alleles A and B respectively computed from average parental Gij or parental genotype counts (i.e. the elements of 
  #               the flj vector).
  #          4-7. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the fjfj matrix (see the 
  #               bottom left of page 4 of Henshall et al. 2014
  #          8. SAMPLE_ID is the individual identifier (all values are zero and represent unknown individuals)
  
  # miss.parent.count: Data frame
  #          PARENT_COMBN_ID     
  #          SNP_ID         
  #          MISS_PARENT_COUNT
  
  #Returns##########################################
  
  # duos.lod: Data frame.
  #          1. OFFSPRING_ID: Offspring identifier from individual identifier in Gij
  #          2. PARENT_COMBN_ID: Combinations of families 
  #          3. SNP_COUNT: The number of SNPS for which none of the sire, dam or offspring genotypes are missing
  #          4. MIN_LOGL: Minimum log-likelihood ratio of SNPs (that were not excluded) in duo.  
  #          5. MAX_LOGL: Maximum log-likelihood ratio of SNPs in duo.
  #          6. RANGE_5_TO_95_LOGL: 5-95 percentile range of log-likelihood ratio
  #          7. LOD: Summed log-likelihood ratios (log odds score) for the offspring, sire, dam combination.
  
  # duos.logl: Data frame.
  #          1. SNP_ID:	SNP identifier
  #          2. OFFSPRING_ID:	Offspring identifier
  #          3. PARENT_COMBN_ID: Combinations of families 
  #          4. LIKE_NUMERATOR:	The likelihood of the offspring, sire and dam duo for the given SNP_  See Equation 3 of Henshall et al. (2014).
  #          5. LIKE_DENOMINATOR:	"The likelihood under the null hypothesis that the offspring is unrelated to the sire and dam, which is constructed by treating both parents as missing" Henshall et al. (2014).
  #          6. LOG_LIKE_RATIO:	Log likelihood ratio of the offspring, sire and dam duo for the given SNP_
  
  print("Running parents.to.pools.lod.fun")
  
  # print("Still running fam.to.pools.lod.fun")
  
  #Loop through pools 
  
  #seed data 
  duos.lod <- NULL
  n.in.pools <- (ncol(parent.combns.by.fam.set.combn)-2)/2
  
  #identify pools 
  pools <- unique(g.d.klj.adj$SAMPLE_ID)
  
  #generatefam.set.combns.by.pool
  fam.set.combns.by.pool <- unique(g.d.klj.adj[,c("SAMPLE_ID", "FAM_SET_COMBN_ID")])
  g.d.pj.adj <- g.d.klj.adj[,colnames(g.d.klj.adj) != "FAM_SET_COMBN_ID"]
  rm(g.d.klj.adj)
  
  #remove duplicated parent combinations within FAM_SET_COMBN_ID
  parent.combns.by.fam.set.combn <- parent.combns.by.fam.set.combn[!duplicated(parent.combns.by.fam.set.combn[,c("FAM_SET_COMBN_ID","PARENT_COMBN_ID")]), ]
  
  #Add MISS_PARENT_SNP_DATA_PROP to miss.parent.count
  if(!is.null(miss.parent.count)) {
    miss.parent.count[,"MISS_PARENT_SNP_DATA_PROP"]  <- miss.parent.count[,"MISS_PARENT_COUNT"] / (n.in.pools * 2)
    miss.parent.count <- miss.parent.count[,c("PARENT_COMBN_ID", "SNP_ID", "MISS_PARENT_SNP_DATA_PROP")]
  }  
  
  for(i in 1:length(pools)) {
    
    pool <- pools[i] #pools identifier    
    print(paste("Computing duo LODs for sample", pool, "-", i, "of", length(pools)))
    
    fam.set.combn <-fam.set.combns.by.pool[fam.set.combns.by.pool[,"SAMPLE_ID"] == pool,"FAM_SET_COMBN_ID"]
    
    #only retain parent combinations for the fam.set.combn
    t.adj.cj <- tclj.adj[tclj.adj[,"PARENT_COMBN_ID"] %in% 
                           parent.combns.by.fam.set.combn[parent.combns.by.fam.set.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn, "PARENT_COMBN_ID"],]
    t.adj.cj <- t.adj.cj[,colnames(t.adj.cj) != "FAM_SET_COMBN_ID"]
    
    nj <- nlj[nlj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    nj <- nj[,colnames(nj) != "FAM_SET_COMBN_ID"]
    
    print("Generating log likelihood numerators")    
    
    # like.num <- merge(t.adj.cj,
    #                   g.d.pj.adj[g.d.pj.adj[,"SAMPLE_ID"] == pool,],
    #                   by = "SNP_ID",
    #                   suffixes = c("_TRANS","_GENO_PROB"),
    #                   all = TRUE
    # )    
    
    like.num <- inner_join(t.adj.cj, 
                           g.d.pj.adj[g.d.pj.adj[,"SAMPLE_ID"] == pool,], 
                           suffix = c("_TRANS","_GENO_PROB"),
                           by = "SNP_ID")
    like.num <- like.num[order(like.num$SNP_ID),]
    
    #Generate list of genotypes
    genotypes <- genotypes.fun(n=(2*n.in.pools))
    
    like.num[,"LIKE_NUMERATOR"] <- rowSums(like.num[,paste(genotypes, "_TRANS", sep="")] * 
                                             like.num[,paste(genotypes, "_GENO_PROB", sep="")])
    
    like.num <- like.num[,c("SNP_ID", "SAMPLE_ID", "PARENT_COMBN_ID", "LIKE_NUMERATOR")]
    
    if (nrow(like.num) > 0) {
      
      print("Generating log likelihood denominators")
      
      #Denominator in likelihood ratio   
      #  like.den <- merge(nj,
      #                    g.d.pj.adj[g.d.pj.adj[,"SAMPLE_ID"] == pool,],
      #                    by = "SNP_ID",
      #                    suffixes = c("_TRANS","_GENO_PROB"),
      #                    all = TRUE
      #  )
      
      like.den <- inner_join(nj, 
                             g.d.pj.adj[g.d.pj.adj[,"SAMPLE_ID"] == pool,], 
                             suffix = c("_TRANS","_GENO_PROB"),
                             by = "SNP_ID")
      
      #Equation 3 of Henshall et al. 2014
      like.den[,"LIKE_DENOMINATOR"] <- rowSums(like.den[,paste(genotypes, "_TRANS", sep="")] * 
                                                 like.den[,paste(genotypes, "_GENO_PROB", sep="")])
      
      like.den <- like.den[,c("SNP_ID", "LIKE_DENOMINATOR")]
      
      print("Generating log likelihood ratios")
      
      #Likelihood ratio
      # tmp.duos <- merge(like.num, like.den, by = "SNP_ID", all.x = TRUE)
      tmp.duos <- left_join(like.num, like.den, by = "SNP_ID")
      tmp.duos$LOG_LIKE_RATIO <- log(tmp.duos$LIKE_NUMERATOR) - log(tmp.duos$LIKE_DENOMINATOR)
      tmp.duos$PARENT_COMBN_ID  <- as.integer(tmp.duos$PARENT_COMBN_ID)
      
      #Add columns MISS_PARENT_SNP_DATA_PROP and MISS_POOL_SNP_DATA_PROP to tmp.duos
      
      if(is.null(miss.parent.count)) {
        tmp.duos[,"MISS_PARENT_SNP_DATA_PROP"] <- NA
      } else {
        #tmp.duos <- merge(tmp.duos, miss.parent.count, by = c("PARENT_COMBN_ID", "SNP_ID"), all.x = TRUE)
        tmp.duos$PARENT_COMBN_ID <- as.integer(tmp.duos$PARENT_COMBN_ID)
        tmp.duos <- left_join(tmp.duos, miss.parent.count,by = c("PARENT_COMBN_ID" , "SNP_ID"))
        tmp.duos <- tmp.duos[order(tmp.duos$PARENT_COMBN_ID),]
      }  
      
      miss.g.d.pj.adj <- g.d.pj.adj[,c("SNP_ID", "SAMPLE_ID")]
      miss.g.d.pj.adj[,"MISS_POOL_SNP_DATA_PROP"] <- as.numeric(is.na(g.d.pj.adj[,3]))
      # tmp.duos <- merge(tmp.duos, miss.g.d.pj.adj, by = c("SNP_ID", "SAMPLE_ID"), all.x = TRUE)
      tmp.duos <- left_join(tmp.duos, miss.g.d.pj.adj, by = c("SNP_ID" , "SAMPLE_ID"))
      tmp.duos <- tmp.duos[order(tmp.duos$SNP_ID),]
      rm(miss.g.d.pj.adj)
      
      #Retain all data in tmp.duos
      duos.logl.all <- tmp.duos
      
      duos.logl <- tmp.duos
      
      print("Generating LODs")
      
      #sum across all snp to get 'log odds score' (LOD)  
      tmp.duos <- aggregate(tmp.duos$LOG_LIKE_RATIO, by = list(tmp.duos$SAMPLE_ID,
                                                               tmp.duos$PARENT_COMBN_ID), 
                            na.rm=T, 
                            FUN = "sum")  
      colnames(tmp.duos) <- c("SAMPLE_ID", "PARENT_COMBN_ID", "LOD")
      #  tmp.duos$SAMPLE_ID <- pool
      tmp.duos <- tmp.duos[,c("SAMPLE_ID", "PARENT_COMBN_ID", "LOD")]
      
      
      #Get maximum LOG_LIKE_RATIO
      
      print("Determining maximum log likelihood ratios")
      max.logl <- suppressWarnings(aggregate(duos.logl$LOG_LIKE_RATIO, by = list(duos.logl$SAMPLE_ID,
                                                                                 duos.logl$PARENT_COMBN_ID), 
                                             na.rm=T, 
                                             FUN = "max"))
      colnames(max.logl) <- c("SAMPLE_ID","PARENT_COMBN_ID", "MAX_LOGL")
      
      #    max.logl <- merge(max.logl, 
      #                      duos.logl[, c("PARENT_COMBN_ID", "LOG_LIKE_RATIO", "SNP_ID")], 
      #                      by.x = c("PARENT_COMBN_ID", "MAX_LOGL"),
      #                      by.y = c("PARENT_COMBN_ID", "LOG_LIKE_RATIO"),
      #                      all.x = TRUE)
      
      tmp <- c("PARENT_COMBN_ID" = "PARENT_COMBN_ID", "MAX_LOGL" = "LOG_LIKE_RATIO")
      max.logl <- left_join(max.logl, 
                            duos.logl[, c("PARENT_COMBN_ID", "LOG_LIKE_RATIO", "SNP_ID")], 
                            by = tmp)
      rm(tmp)
      
      colnames(max.logl) <- c("SAMPLE_ID", "PARENT_COMBN_ID", "MAX_LOGL", "MAX_LOGL_SNP")
      #    max.logl$SAMPLE_ID   <- pool
      #max.logl           <- max.logl[,c("SAMPLE_ID", "PARENT_COMBN_ID", "MAX_LOGL", "MAX_LOGL_SNP")]
      
      #Get minimum LOG_LIKE_RATIO
      
      print("Determining minimum log likelihood ratios")
      
      min.logl <- suppressWarnings(aggregate(duos.logl$LOG_LIKE_RATIO, by = list(duos.logl$SAMPLE_ID,
                                                                                 duos.logl$PARENT_COMBN_ID), 
                                             na.rm=T, 
                                             FUN = "min")) 
      colnames(min.logl) <- c("SAMPLE_ID","PARENT_COMBN_ID", "MIN_LOGL")
      
      #  min.logl <- merge(min.logl, 
      #                    duos.logl[, c("PARENT_COMBN_ID", "LOG_LIKE_RATIO", "SNP_ID")], 
      #                    by.x = c("PARENT_COMBN_ID", "MIN_LOGL"),
      #                    by.y = c("PARENT_COMBN_ID", "LOG_LIKE_RATIO"),
      #                    all.x = TRUE)
      
      tmp <- c("PARENT_COMBN_ID" = "PARENT_COMBN_ID", "MIN_LOGL" = "LOG_LIKE_RATIO")
      min.logl <- left_join(min.logl, 
                            duos.logl[, c("PARENT_COMBN_ID", "LOG_LIKE_RATIO", "SNP_ID")], 
                            by = tmp)
      rm(tmp)
      
      colnames(min.logl) <- c("SAMPLE_ID","PARENT_COMBN_ID", "MIN_LOGL",  "MIN_LOGL_SNP")
      # min.logl$SAMPLE_ID <- pool
      #min.logl <- min.logl[,c("SAMPLE_ID", "PARENT_COMBN_ID", "MIN_LOGL", "MIN_LOGL_SNP")]
      
      #Get 5 to 95 quantile range LOG_LIKE_RATIO
      
      print("Determining 5 to 95 quantile range of log likelihood ratios")
      
      ipr.fun <- function(x) {
        quantile(x, c(0.95), na.rm = TRUE) - quantile(x, c(0.05), na.rm = TRUE) 
      }
      iqr.logl <- aggregate(duos.logl$LOG_LIKE_RATIO, by = list(duos.logl$SAMPLE_ID,
                                                                duos.logl$PARENT_COMBN_ID), 
                            #    na.rm=T, 
                            FUN = ipr.fun)  
      colnames(iqr.logl) <- c("SAMPLE_ID", "PARENT_COMBN_ID", "RANGE_5_TO_95_LOGL")
      #  iqr.logl$SAMPLE_ID <- pool
      #  iqr.logl <- iqr.logl[,c("SAMPLE_ID", "PARENT_COMBN_ID", "RANGE_5_TO_95_LOGL")]
      
      print("Determining missing log likelihood ratios")
      
      #geg columns MISS_PARENT_SNP_DATA_PROP and MISS_POOL_SNP_DATA_PROP and NO_MISS_PARENT_OR_POOL_PROP for merge with tmp.duos
      miss.logl <- aggregate(cbind(duos.logl[,c("MISS_PARENT_SNP_DATA_PROP", "MISS_POOL_SNP_DATA_PROP")],
                                   duos.logl[,"MISS_PARENT_SNP_DATA_PROP"] == 0 & duos.logl[,"MISS_POOL_SNP_DATA_PROP"] == 0), 
                             by = list(duos.logl[,c("PARENT_COMBN_ID")], duos.logl[,c("SAMPLE_ID")]),
                             FUN="mean")
      colnames(miss.logl) <- c("PARENT_COMBN_ID", "SAMPLE_ID", "MISS_PARENT_SNP_DATA_PROP", "MISS_POOL_SNP_DATA_PROP", 
                               "NO_MISS_PARENT_OR_POOL_PROP")
      
      #Merge
      tmp.duos <- merge(iqr.logl, tmp.duos, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.y = TRUE)
      tmp.duos <- merge(max.logl, tmp.duos, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.y = TRUE)
      tmp.duos <- merge(min.logl, tmp.duos, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.y = TRUE)
      tmp.duos <- merge(miss.logl, tmp.duos, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.y = TRUE)      
      
      duos.lod <- rbind(duos.lod,tmp.duos)
      #     rm(like.num, like.den, tmp.duos)
    } #END if (nrow(like.num) > 0) {
  } #END   for(i in 1:length(pools)) {
  
  return(list(duos.lod = duos.lod, duos.logl = duos.logl.all))
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

most.like.parents.duo.fun <- function(duos.lod,
                                      parent.combns) {
  
  #Description: Identifies the family that each inidividual is most likely to belong to.
  
  #  Args:
  # duos.lod: Data frame.
  #          1. FAM_COMBN_ID: Family combination identifier from parent.combns.by.fam.set.combn
  #          2. SAMPLE_ID: Sample identifier
  
  #          ????. SNP_COUNT: The number of SNPS for which none of the sire, dam or offspring genotypes are missing
  
  #          3. MIN_LOGL: Minimum log-likelihood ratio of SNPs in duo
  #          4. MIN_LOGL_SNP: SNP associated with MIN_LOGL
  #          5. MAX_LOGL: Maximum log-likelihood ratio of SNPs in duo
  #          6. MAX_LOGL_SNP: SNP associated with MAX_LOGL
  #          7. RANGE_5_TO_95_LOGL: RANGE_5_TO_95_LOGL: 5-95 percentile range of log-likelihood ratio
  #          8. LOD: Summed log-likelihood ratios (log odds score) for the FAM_COMBN_ID and SAMPLE_ID duo.
  #          9. FAMILY_ID_1: Family 1 in FAM_COMBN_ID in parent.combns.by.fam.set.combn
  #          etc
  
  # parent.combns.by.fam.set.combn
  
  
  
  
  #  Returns:
  #parent.assign: Data frame.
  #          1. FAM_COMBN_ID: Family combination identifier from parent.combns.by.fam.set.combn with the highest LOD in combination with SAMPLE_ID
  #          2. SAMPLE_ID: Sample identifier
  
  #          ????. SNP_COUNT: The number of SNPS for which none of the sire, dam or offspring genotypes are missing
  
  #          3. MIN_LOGL: Minimum log-likelihood ratio of SNPs in duo
  #          4. MIN_LOGL_SNP: SNP associated with MIN_LOGL
  #          5. MAX_LOGL: Maximum log-likelihood ratio of SNPs in duo
  #          6. MAX_LOGL_SNP: SNP associated with MAX_LOGL
  #          7. RANGE_5_TO_95_LOGL: RANGE_5_TO_95_LOGL: 5-95 percentile range of log-likelihood ratio
  #          8. LOD: Summed log-likelihood ratios (log odds score) for the FAM_COMBN_ID and SAMPLE_ID duo.
  #          9. FAMILY_ID_1: Family 1 in FAM_COMBN_ID in parent.combns.by.fam.set.combn
  #          10 FAM_1_DELTA_LOD: Highest LOD of FAM_COMBN SAMPLE_ID if FAMILY_ID_1 excluded
  #          etc
  
  print("Running most.like.parents.duo.fun")
  
  #remove FAM_COMBN_ID containing duplicated PARENT_COMBN_ID
  tmp <- parent.combns[duplicated(parent.combns[,"PARENT_COMBN_ID"]),
                       "PARENT_COMBN_ID"]
  
  if(length(tmp) != 0) {
    parent.combns[parent.combns[,"PARENT_COMBN_ID"] %in% tmp,grep("FAM",colnames(parent.combns))] <- NA
  }
  parent.combns.unambiguous <- parent.combns
  #remove PARENT_IDs
  parent.combns.unambiguous <- parent.combns.unambiguous[,colnames(parent.combns.unambiguous)[c(grep("FAM", colnames(parent.combns.unambiguous)), 
                                                                                                grep("PARENT_COMBN_ID", colnames(parent.combns.unambiguous)))]
                                                         ]
  
  #Obtain most likely family (i.e. maximum LOD) for each offspring
  most.like           <- aggregate(duos.lod$LOD, by = list(duos.lod$SAMPLE_ID), na.rm=T, FUN = "max")      
  colnames(most.like) <- c("SAMPLE_ID", "LOD")
  
  # most.like <- merge(most.like, duos.lod, by = c("SAMPLE_ID", "LOD"), all.x = TRUE)
  most.like <- left_join(most.like, duos.lod, by = c("SAMPLE_ID", "LOD"))
  # most.like <- merge(most.like, parent.combns.unambiguous, by = c("PARENT_COMBN_ID"), all.x = TRUE)
  most.like <- left_join(most.like, parent.combns.unambiguous, by = c("PARENT_COMBN_ID"))
  most.like <- unique(most.like)
  
  #Remove most likely duo to obtain second most likely
  most.like.parents.list <- most.like[,c("SAMPLE_ID", "PARENT_COMBN_ID")]
  # most.like.parents.list <- merge(most.like.parents.list, 
  #                                  parent.combns[,colnames(parent.combns)[grep("PARENT", colnames(parent.combns))] ]  , 
  #                                  by = c("PARENT_COMBN_ID"), all.x = TRUE)
  most.like.parents.list <- left_join(most.like.parents.list, 
                                      parent.combns[,colnames(parent.combns)[grep("PARENT", colnames(parent.combns))] ], 
                                      by = c("PARENT_COMBN_ID"))
  #cycle for n.in.pools
  n.in.pools <- (ncol(parent.combns[,colnames(parent.combns)[grep("PARENT", colnames(parent.combns))] ] )-1)/2
  
  #identify samples with ambiguous parent allocation
  most.like.parents.list <- unique(most.like.parents.list)
  ambig.assign.list <- NULL
  ambig.assign.list <- unique(most.like.parents.list[duplicated(most.like.parents.list$SAMPLE_ID),"SAMPLE_ID"])
  
  #remove duplicate sample ids and then overwrite outputs for samples in ambig.assign.list at end
  most.like.parents.list <- most.like.parents.list[!duplicated(most.like.parents.list$SAMPLE_ID),]
  
  for(i in 1:(2*n.in.pools)) {
    #identify parent i for each sample
    
    tmp.remove           <- most.like.parents.list[,c("SAMPLE_ID",paste("PARENT_",i,sep = ""))]
    
    colnames(tmp.remove) <- c("SAMPLE_ID", "PARENT_REMOVE")
    tmp.remove[,"PARENT_REMOVE_COUNT"] <- rowSums(most.like.parents.list[,paste("PARENT_",1:(2*n.in.pools),sep = "")] == 
                                                    tmp.remove[,"PARENT_REMOVE"])
    tmp.remove <- unique(tmp.remove)
    
    # duos.lod.removed     <- merge(duos.lod, tmp.remove, by = "SAMPLE_ID", all.x = TRUE)
    duos.lod.removed <- left_join(duos.lod, tmp.remove, by = "SAMPLE_ID")
    
    duos.lod.removed[,"REMOVE"] <- 
      rowSums(duos.lod.removed[,paste("PARENT_",1:(2*n.in.pools),sep = "")] == 
                duos.lod.removed[,"PARENT_REMOVE"]) >= 
      duos.lod.removed[,"PARENT_REMOVE_COUNT"]
    
    duos.lod.removed <- duos.lod.removed[!duos.lod.removed$REMOVE,]
    
    if(nrow(duos.lod.removed) > 0) { 
      #Obtain second most likely family (i.e. maximum LOD) for each offspring
      second.most.like <-  aggregate(duos.lod.removed$LOD, 
                                     by = list(duos.lod.removed$SAMPLE_ID), 
                                     na.rm=T, FUN = "max") 
      colnames(second.most.like) <- c("SAMPLE_ID", paste("PARENT_",i, "_ALT_LOD",sep = ""))
      second.most.like <- second.most.like[!duplicated(second.most.like),]
      # second.most.like <- merge(tmp.remove, second.most.like, by = "SAMPLE_ID", all.x = TRUE) #in case not all samples in second.most.like
      second.most.like <- left_join(tmp.remove, second.most.like, by = "SAMPLE_ID")
      
      
      #Identify alternative parent
      alt.parents <- merge(second.most.like,
                           duos.lod.removed, 
                           by.x = c("SAMPLE_ID", paste("PARENT_",i, "_ALT_LOD",sep = "")),
                           by.y  = c("SAMPLE_ID", "LOD"), 
                           all.x = TRUE)[,c("SAMPLE_ID", "PARENT_COMBN_ID", paste("PARENT_",1:(2*n.in.pools),sep = ""))]
      colnames(alt.parents) <- c("SAMPLE_ID", "PARENT_COMBN_ID", paste("ALT_PARENT_",1:(2*n.in.pools),sep = ""))
      
      multi.alt.parents <-  aggregate(!is.na(alt.parents$SAMPLE_ID), 
                                      by = list(alt.parents$SAMPLE_ID), 
                                      na.rm=T, FUN = "sum") 
      alt.parents <- alt.parents[!duplicated(alt.parents[,"SAMPLE_ID"]),]
      alt.parent <- NULL
      
      for(r in 1:nrow(alt.parents)) {
        
        tmp.1 <- aggregate(t(!is.na(most.like.parents.list[r,3:(2*n.in.pools+2)])), 
                           by = list(as.integer(most.like.parents.list[r,3:(2*n.in.pools+2)])), 
                           na.rm=T, FUN = "sum") 
        
        tmp.2 <-  aggregate(t(!is.na(alt.parents[r,3:(2*n.in.pools+2)])), 
                            by = list(as.integer(alt.parents[r,3:(2*n.in.pools+2)])), 
                            na.rm=T, FUN = "sum") 
        if (nrow(tmp.2) > 0) {
          tmp.3 <- merge(tmp.1, tmp.2, by = 1, all = TRUE)
        } else {
          tmp.3 <- tmp.1
          tmp.3[,3] <- NA
        }
        
        tmp.3[is.na(tmp.3[,1]),1] <- 0
        tmp.3[is.na(tmp.3[,2]),2] <- 0
        tmp.3[is.na(tmp.3[,3]),3] <- 0
        
        tmp <- tmp.3[tmp.3[,3]-tmp.3[,2] == 1,1]
        if(length(tmp) != 1) {tmp <- NA} #more than one parent changed or no alternatives
        rm(tmp.1, tmp.2, tmp.3)
        # tmp <- setdiff(alt.parents[r,3:(2*n.in.pools+2)], most.like.parents.list[r,3:(2*n.in.pools+2)])
        #tmp <- c(alt.parents[r,"SAMPLE_ID"], tmp)
        alt.parent <- c(alt.parent,tmp)
      }
      
      alt.parent <-  unlist(alt.parent)
      alt.parent <-  cbind(multi.alt.parents, alt.parent)
      colnames(alt.parent) <- c("SAMPLE_ID", "COUNT", paste("ALT_PARENT_", i, sep=""))
      alt.parent[alt.parent[,"COUNT"] > 1,paste("ALT_PARENT_", i, sep="")] <- NA
      alt.parent <- alt.parent[,c("SAMPLE_ID", paste("ALT_PARENT_", i, sep=""))]
      
      # most.like <- merge(most.like, alt.parent, by = "SAMPLE_ID", all.x = TRUE)
      most.like <- left_join(most.like, alt.parent, by = "SAMPLE_ID")
      
      alt.combn <- alt.parents[,c("SAMPLE_ID", "PARENT_COMBN_ID")]
      colnames(alt.combn) <- c("SAMPLE_ID", paste("ALT_PARENT_COMBN_", i, sep=""))
      # most.like <- merge(most.like, alt.combn, by = "SAMPLE_ID", all.x = TRUE)
      most.like <- left_join(most.like, alt.combn, by = "SAMPLE_ID")
      
      rm(tmp, alt.parents,multi.alt.parents, alt.parent, alt.combn)
      
      #Get delta LOD
      second.most.like <- second.most.like[,c("SAMPLE_ID", paste("PARENT_",i, "_ALT_LOD",sep = ""))]
      #most.like <- merge(most.like, second.most.like, by = "SAMPLE_ID", all.x = TRUE)
      most.like <- left_join(most.like, second.most.like, by = "SAMPLE_ID")
      
      most.like[,paste("PARENT_",i, "_DELTA_LOD",sep = "")] <-
        most.like[,"LOD"] - most.like[,paste("PARENT_",i, "_ALT_LOD",sep = "")]
      most.like <- most.like[,colnames(most.like) !=  paste("PARENT_",i, "_ALT_LOD",sep = "")]
      
    } else {
      most.like[,paste("PARENT_",i, "_DELTA_LOD",sep = "")] <- NA
      most.like[,paste("ALT_PARENT_",i, sep = "")] <- NA
      most.like[,paste("ALT_PARENT_COMBN_", i, sep = "")] <- NA
    }
  }
  
  #column order
  
  #Alternate family ids
  parent.combns.unambiguous <- parent.combns.unambiguous[,c("PARENT_COMBN_ID", "FAM_COMBN_ID")]
  colnames(parent.combns.unambiguous) <- c("PARENT_COMBN_ID", "ALT_FAM_COMBN_ID")
  parent.combns.unambiguous <- parent.combns.unambiguous[!is.na(parent.combns.unambiguous[,"ALT_FAM_COMBN_ID"]),]
  parent.combns.unambiguous <- unique(parent.combns.unambiguous)
  
  for(i in 1:(2*n.in.pools)) {
    
    most.like <- merge(most.like, 
                       parent.combns.unambiguous, 
                       by.x = paste("ALT_PARENT_COMBN_",i, sep=""),
                       by.y = "PARENT_COMBN_ID",
                       all.x = TRUE)
    colnames(most.like)[colnames(most.like) == "ALT_FAM_COMBN_ID"] <-  paste("ALT_FAM_COMBN_",i, sep="")
  }
  
  most.like <- most.like[,c(colnames(duos.lod), 
                            "FAM_COMBN_ID",
                            paste("FAMILY_ID_",1:n.in.pools,sep = ""),
                            paste("PARENT_",1:(2*n.in.pools), "_DELTA_LOD",sep = ""),
                            paste("ALT_PARENT_",1:(2*n.in.pools),sep = ""),
                            paste("ALT_PARENT_COMBN_",1:(2*n.in.pools),sep = ""),
                            paste("ALT_FAM_COMBN_",1:(2*n.in.pools),sep = ""))
                         ]
  
  #Remove data for samples in ambig.assign.list
  most.like <- most.like[!most.like[,"SAMPLE_ID"] %in% ambig.assign.list,]
  #Add row of NA for samples in ambig.assign.list
  # most.like <- merge(most.like, data.frame(SAMPLE_ID = ambig.assign.list), by = "SAMPLE_ID", all = TRUE)
  tmp <- data.frame(SAMPLE_ID = ambig.assign.list)
  most.like <- left_join(most.like, tmp, by = "SAMPLE_ID")
  
  most.like <- most.like[order(most.like[,"SAMPLE_ID"]),]
  
  return(most.like)
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

plot.lod.vs.LogL.range.fun <- function(lod, meth) {
  
  #generates scatterplot
  
  #Args##########################################
  
  # trios.lod: Data frame.  Output of trios.lod.fun
  #          1. OFFSPRING_ID: Offspring identifier from individual identifier in Gij
  #          9. RANGE_5_TO_95_LOGL: 5-95 percentile range of log-likelihood ratio
  #          10. LOD: Summed log-likelihood ratios (log odds score) for the offspring, sire, dam combination.
  
  #meth = current method
  
  #Returns##########################################
  
  # scatterplot:  Plot of LOD vs RANGE_5_TO_95_LOGL returned 
  #               in a directory named "lod.scatter" in the working directory. 
  
  print("Running plot.lod.vs.LogL.range.fun")
  
  #Check that all headings are present in inputs  
  #  if(sum(c("OFFSPRING_ID",	"SIRE_ID",	"DAM_ID",	"FAMILY_ID",	"SNP_COUNT",	"SNP_EXCLUDED",	"MIN_LOGL",	"MIN_LOGL_SNP",
  #           "MAX_LOGL",	"MAX_LOGL_SNP",	"RANGE_5_TO_95_LOGL",	"LOD") %in% colnames(trios.lod)) != 12) {
  #    stop("indiv.snp.dat.indiv input must be a data frame containing the following headings: OFFSPRING_ID, SIRE_ID, DAM_ID, FAMILY_ID, SNP_COUNT, SNP_EXCLUDED, MIN_LOGL, MIN_LOGL_SNP, MAX_LOGL, MAX_LOGL_SNP, RANGE_5_TO_95_LOGL, LOD")
  #  }
  
  #Name columns and assign class
  
  #  trios.lod$OFFSPRING_ID  <- as.character(trios.lod$OFFSPRING_ID)
  #  trios.lod$RANGE_5_TO_95_LOGL    <- as.numeric(trios.lod$RANGE_5_TO_95_LOGL)
  #  trios.lod$LOD    <- as.numeric(trios.lod$LOD)
  
  wd <- getwd()
  #Generate directory for plots if does not exist
  dir.create(file.path(wd, "lod.scatter"), showWarnings = FALSE)
  setwd(file.path(wd, "lod.scatter"))
  
  #Loop through offspring
  for(sample in unique(lod[,1])) {
    
    tmp.lod <- lod[lod[,1] == sample,]
    
    #Define plot variables and data
    tmp.lod$X.AXIS <- tmp.lod$LOD
    x.title <- "Log odds (LOD) score"
    tmp.lod$Y.AXIS <- tmp.lod$RANGE_5_TO_95_LOGL    
    y.title <- "5-95 percentile range of log-likelihood ratio"
    main.title <- sample
    
    png(filename = paste(sample, meth, ".png", sep=""))
    
    #Generate plot
    plot(tmp.lod$X.AXIS,tmp.lod$Y.AXIS,
         type="p",pch = 16, cex = 0.5,
         main= main.title,
         xlab = x.title,
         ylab = y.title)
    
    dev.off()
  }
  setwd(wd)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

prelim.ml.discrete.assigned.genos.fun <- function(method,
  snp.dat.indiv,
  snp.dat.pools,
  parent.combns.by.fam.set.combn,
  flj,
  snp.error,
  nlj,
  fam.set.combns.by.pool) {
  
  print("Running prelim.ml.discrete.assigned.genos.fun")
  
  Dij <- Dij.from.snp.dat.fun(snp.dat.indiv = snp.dat.indiv) 
  
  tclj <- NULL
  for(fam.set.combn in unique(flj$FAM_SET_COMBN_ID)) {
    
    tmp.fj <- flj[flj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    tmp.parent.combns <- parent.combns.by.fam.set.combn[parent.combns.by.fam.set.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn,] 
    
    tmp.parents <- tmp.parent.combns[,!colnames(tmp.parent.combns) %in% c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID")]
    tmp.parents <- unique(as.vector(as.matrix(tmp.parents)))
    
    tmp.tij <- Dij[Dij[,"SAMPLE_ID"] %in% tmp.parents, ]
    tmp.tclj <- tcj.fun(tij = tmp.tij, 
                        parent.combns = tmp.parent.combns, 
                        fj = tmp.fj)
    tmp.tclj$FAM_SET_COMBN_ID <- fam.set.combn
    tclj <- rbind(tclj, tmp.tclj)
  }
  rm(tmp.fj, tmp.parent.combns, tmp.parents, tmp.tij, tmp.tclj)
  
  miss.parent.count <- tclj[,c("PARENT_COMBN_ID", "SNP_ID", "MISS_PARENT_COUNT")]
  tclj               <- tclj[,colnames(tclj) != "MISS_PARENT_COUNT"] 
  
  if(!"B_ALLELE" %in% colnames(snp.dat.pools)) {
  snp.dat.pools <- left_join(snp.dat.pools, 
                             unique(snp.dat.indiv[snp.dat.indiv[,"SNP_ID"] %in% snp.dat.pools[,"SNP_ID"],
                                                          c("SNP_ID", "A_ALLELE", "B_ALLELE")]), 
                             by = "SNP_ID")
  }
  
  dkj <- dkj.from.snp.dat.pools.fun(snp.dat.pools = snp.dat.pools) 
  
  if("Discrete" %in% method) {
  g.d.klj.adj <- adj.geno.prob.fun(lhs = merge(dkj,fam.set.combns.by.pool, by = "SAMPLE_ID", all = TRUE),
                                   snp.error=snp.error,
                                   nlj = nlj)
  } else {
    g.d.klj.adj <- NULL
  }
  
  return(list(  Dij = Dij, 
                tclj = tclj,
                miss.parent.count = miss.parent.count,
                dkj = dkj,
                g.d.klj.adj = g.d.klj.adj))
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

Dij.from.snp.dat.fun <- function(snp.dat.indiv) {
  
  #Returns########################################## the Dij matrix for each individual and SNP according to the descrition on page 5 of Henshall et al. 2014.   
  #Discrete genotypes assigned (1 for homozygote, 0.5 for heterozygote)
  #to an element of Dij based on GENOTYPE in snp.dat.indiv.
  
  #Args##########################################
  # snp.dat.indiv: Data frame
  #              SAMPLE_ID is the individual identifier
  #              SNP_ID   is the SNP identifier
  #              A_ALLELE  is the base represented by allele A
  #              B_ALLELE  is the base represented by allele B
  #              GENOTYPE is the SNP genotype call
  
  #Returns##########################################
  # D.matrices:     Data frame 
  #              1.   SAMPLE_ID is the individual identifier 
  #              2.   SNP_ID is the SNP identifier
  #              3-6. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the Dij matrix (see the 
  #                   top left of page 5 of Henshall et al. 2014
  #              7-8. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #                   for alleles A and B respectively computed from Dij (i.e. the elements of 
  #                   the transmission (Tij) vector, Equation 2 of Henshall et al. 2014).
  
  print("Running Dij.from.snp.dat.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID","A_ALLELE", "B_ALLELE", "GENOTYPE") %in% colnames(snp.dat.indiv)) != 5) {
    stop("snp.dat.indiv input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, A_ALLELE, B_ALLELE, GENOTYPE")
  }
  
  #Name columns and assign class
  #colnames(snp.dat.indiv) <- c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B", "A_ALLELE", "B_ALLELE", "GENOTYPE")
  snp.dat.indiv$SAMPLE_ID  <- as.integer(snp.dat.indiv$SAMPLE_ID)
  snp.dat.indiv$SNP_ID    <- as.character(snp.dat.indiv$SNP_ID)
  snp.dat.indiv$A_ALLELE   <- as.character(snp.dat.indiv$A_ALLELE)
  snp.dat.indiv$B_ALLELE   <- as.character(snp.dat.indiv$B_ALLELE)
  snp.dat.indiv$GENOTYPE  <- as.character(snp.dat.indiv$GENOTYPE)
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "A_ALLELE", "B_ALLELE", "GENOTYPE")]
  
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
            paste(snp.dat.indiv[,"B_ALLELE"], snp.dat.indiv[,"A_ALLELE"], sep = ""),"ALLELES"] <- "BA"
  
  #Change missing values back to NA
  snp.dat.indiv[snp.dat.indiv[,"GENOTYPE"] == "No.call", "GENOTYPE"] <- NA
  
  snp.dat.indiv[!is.na(snp.dat.indiv[,"ALLELES"]) , "AA_GENO_PROB"] <- 0
  snp.dat.indiv[!is.na(snp.dat.indiv[,"ALLELES"]) , "AB_GENO_PROB"] <- 0
  snp.dat.indiv[!is.na(snp.dat.indiv[,"ALLELES"]) , "BA_GENO_PROB"] <- 0
  snp.dat.indiv[!is.na(snp.dat.indiv[,"ALLELES"]) , "BB_GENO_PROB"] <- 0
  
  #Assign descrete genotypes.  Note that it is assumed that Gij and Dij are symetrical (i.e. AB = BA = 0.5). 
  # See Henshall et al. 2014 in the top left of page 5.
  
  snp.dat.indiv[snp.dat.indiv[,"ALLELES"] == "AA" & !is.na(snp.dat.indiv[,"ALLELES"]) , "AA_GENO_PROB"] <- 1
  
  snp.dat.indiv[snp.dat.indiv[,"ALLELES"] == "AB" & !is.na(snp.dat.indiv[,"ALLELES"]) , "AB_GENO_PROB"] <- 1/2
  snp.dat.indiv[snp.dat.indiv[,"ALLELES"] == "BA" & !is.na(snp.dat.indiv[,"ALLELES"]) , "AB_GENO_PROB"] <- 1/2
  
  snp.dat.indiv["BA_GENO_PROB"] <- snp.dat.indiv["AB_GENO_PROB"] 
  
  snp.dat.indiv[snp.dat.indiv[,"ALLELES"] == "BB" & !is.na(snp.dat.indiv[,"ALLELES"]) , "BB_GENO_PROB"] <- 1
  
  #Equation 2 of Henshall et al. 2014 - T vector
  snp.dat.indiv$A_TRANS_PROB <- (snp.dat.indiv$AA_GENO_PROB * 2 + snp.dat.indiv$AB_GENO_PROB + snp.dat.indiv$BA_GENO_PROB)/2  #(sum row 1 + sum col 1) / 2
  snp.dat.indiv$B_TRANS_PROB <- (snp.dat.indiv$BB_GENO_PROB * 2 + snp.dat.indiv$AB_GENO_PROB + snp.dat.indiv$BA_GENO_PROB)/2  #(sum row 2 + sum col 2) / 2
  
  #retain only relevant columns
  snp.dat.indiv <- snp.dat.indiv[,c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB", "A_TRANS_PROB", "B_TRANS_PROB")]
  
  return(snp.dat.indiv)
  
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

dkj.from.snp.dat.pools.fun <- function(snp.dat.pools) {
  
  #Returns########################################## the Dij matrix for each individual and SNP according to the descrition on page 5 of Henshall et al. 2014.   
  #Discrete genotypes assigned (1 for homozygote, 0.5 for heterozygote)
  #to an element of Dij based on GENOTYPE in snp.dat.pools.
  
  #Args##########################################
  # snp.dat.pools: Data frame
  #              SAMPLE_ID is the individual identifier
  #              SNP_ID   is the SNP identifier
  #              A_ALLELE  is the base represented by allele A
  #              B_ALLELE  is the base represented by allele B
  #              GENOTYPE is the SNP genotype call
  
  #Returns##########################################
  # D.matrices:     Data frame 
  #              1.   SAMPLE_ID is the individual identifier 
  #              2.   SNP_ID is the SNP identifier
  
  print("Running dkj.from.snp.dat.pools.fun")
  
  n.in.pools <- max(nchar(snp.dat.pools[,"GENOTYPE"]), na.rm = TRUE)/2
  
  #Change missing values to string of zeros
  tmp.snp.dat.pools <- snp.dat.pools
  tmp.snp.dat.pools[is.na(tmp.snp.dat.pools[,"GENOTYPE"]), "GENOTYPE"] <- paste(as.character(rep(0,n.in.pools*2)),collapse = "")
  
  #count of B allele
  tmp <- t(matrix(unlist(strsplit(tmp.snp.dat.pools[,"GENOTYPE"],"")), 
                  ncol = nrow(tmp.snp.dat.pools)))
  count.b.allele <- rowSums(tmp == tmp.snp.dat.pools[,"B_ALLELE"], na.rm = TRUE)
  count.b.allele[tmp[,1] == 0] <- NA
  
  #generate lambda.kj - unsorted genotype probability vector
  geno.prob <- matrix(0, ncol = (n.in.pools*2+1), nrow = nrow(snp.dat.pools))
  for (i in 1:nrow(geno.prob)) {
    geno.prob[i,(count.b.allele[i]+1)] <- 1
  }
  
  geno.prob[rowSums(geno.prob) == 0,] <- NA
  
  colnames(geno.prob) <- paste(genotypes.fun(n.in.pools*2), "_LAMBDA", sep = "")
  
  lambda.kj <- cbind(snp.dat.pools[,c("SAMPLE_ID", "SNP_ID")], geno.prob)
  
  lambda.kj$ALLELIC_PROP_POOL <- count.b.allele/(n.in.pools*2)
  
  rho.inv <- rho.inv.fun(n.in.pools)
  
  print("Still running dkj.from.snp.dat.pools.fun")
  
  dkj <- lambda.kj
  colnames(dkj) <- gsub("_LAMBDA", "", colnames(dkj))
  dkj <- dkj[,-ncol(dkj)]
  
  for(geno in rho.inv[,"GENOTYPE"]) {
    dkj[,geno] <- dkj[,geno] * rho.inv[rho.inv[,"GENOTYPE"]==geno,"RHO_INV"]
  }
  
  return(dkj)
  
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

prelim.ml.discrete.geno.probs.fun <- function(
  Gij,
  snp.dat.indiv,
  snp.dat.pools,
  threshold.indiv,
  threshold.pools,
  parent.combns.by.fam.set.combn,
  flj,
  gkj,
  lambda.kj,
  snp.error,
  nlj,
  fam.set.combns.by.pool) {
  
  Dij <- Dij.from.Gij.fun(Gij = Gij, threshold.indiv = threshold.indiv) 
  
  tclj <- NULL
  for(fam.set.combn in unique(flj$FAM_SET_COMBN_ID)) {
    
    tmp.fj <- flj[flj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    tmp.parent.combns <- parent.combns.by.fam.set.combn[parent.combns.by.fam.set.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn,] 
    
    tmp.parents <- tmp.parent.combns[,!colnames(tmp.parent.combns) %in% c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID")]
    tmp.parents <- unique(as.vector(as.matrix(tmp.parents)))
    
    tmp.tij <- Dij[Dij[,"SAMPLE_ID"] %in% tmp.parents, ]
    tmp.tclj <- tcj.fun(tij = tmp.tij, 
                        parent.combns = tmp.parent.combns, 
                        fj = tmp.fj)
    tmp.tclj$FAM_SET_COMBN_ID <- fam.set.combn
    tclj <- rbind(tclj, tmp.tclj)
  }
  rm(tmp.fj, tmp.parent.combns, tmp.parents, tmp.tij, tmp.tclj)
  
  miss.parent.count <- tclj[,c("PARENT_COMBN_ID", "SNP_ID", "MISS_PARENT_COUNT")]
  tclj               <- tclj[,colnames(tclj) != "MISS_PARENT_COUNT"] 
  
  dkj <- dkj.from.gkj.fun(lambda.kj = lambda.kj, threshold.pools = threshold.pools) 
  
  g.d.klj.adj <- adj.geno.prob.fun(lhs= merge(dkj,fam.set.combns.by.pool, by = "SAMPLE_ID", all = TRUE),
                                   snp.error=snp.error,
                                   nlj = nlj)
  
  return(list(Dij = Dij,
              tclj = tclj,
              miss.parent.count = miss.parent.count,
              dkj = dkj,
              g.d.klj.adj = g.d.klj.adj))
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

Dij.from.Gij.fun <- function(Gij, threshold.indiv) {
  
  #Description: Returns the Dij matrix for each individual and SNP according to the description on page 5 of 
  #Henshall et al. 2014 from the Gij matrices.  Discrete genotypes assigned (1 for homozygote, 0.5 for heterozygote)
  #to an element of Dij if the corresponding element of Gij greater than a threshold.indiv.
  
  #Args##########################################
  
  #Gij: Data frame.  Output of Gij.fun 
  #          1.   SAMPLE_ID is the individual identifier (assumed to be offspring if not 
  #               listed as a sire or dam in 'fams' input)
  #          2.   SNP_ID is the SNP identifier
  #          3-6. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the Gij  matrix (see the 
  #               top left of page 4 of Henshall et al. 2014
  
  #threshold.indiv: Number.  Genotype called if AA_GENO_PROB, AB_GENO_PROB + BA_GENO_PROB or BB_GENO_PROB greater than this value
  
  #Returns##########################################
  # D.matrices:     Data frame 
  #              1.   SAMPLE_ID is the individual identifier 
  #              2.   SNP_ID is the SNP identifier
  #              3-6. AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB and BB_GENO_PROB: elements of the Dij matrix (see the 
  #                   top left of page 5 of Henshall et al. 2014
  #              7-8. A_TRANS_PROB, B_TRANS_PROB are the probabilities of allele transmission 
  #                   for alleles A and B respectively computed from Dij (i.e. the elements of 
  #                   the transmission (Tij) vector, Equation 2 of Henshall et al. 2014).
  
  print("Running Dij.from.Gij.fun")
  
  #Check that all headings are present in inputs  
  if(sum(c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", 
           "BB_GENO_PROB") %in% 
         colnames(Gij)) != 6) {
    stop("Gij input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, AA_GENO_PROB, AB_GENO_PROB, BA_GENO_PROB, BB_GENO_PROB")
  }
  
  #Name columns and assign class
  #colnames(Gij)     <- c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")
  Gij$SAMPLE_ID     <- as.integer(Gij$SAMPLE_ID)
  Gij$SNP_ID       <- as.character(Gij$SNP_ID)
  Gij$AA_GENO_PROB <- as.numeric(Gij$AA_GENO_PROB)
  Gij$AB_GENO_PROB <- as.numeric(Gij$AB_GENO_PROB)
  Gij$BA_GENO_PROB <- as.numeric(Gij$BA_GENO_PROB)
  Gij$BB_GENO_PROB <- as.numeric(Gij$BB_GENO_PROB)
  Gij <- Gij[,c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")]
  
  #Check that probabilities between 0 and 1
  if(
    sum(
      Gij$AA_GENO_PROB > 1 |
      Gij$AB_GENO_PROB > 1 |
      Gij$BA_GENO_PROB > 1 |
      Gij$BB_GENO_PROB > 1 |
      
      Gij$AA_GENO_PROB < 0 |
      Gij$AB_GENO_PROB < 0 |
      Gij$BA_GENO_PROB < 0 |
      Gij$BB_GENO_PROB < 0 
      , na.rm = T) != 0
  ) {
    stop("Probabilities in Gij must be between 0 and 1 inclusive")
  }
  
  #Check that "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB" sum to one
  if(
    sum(
      round((Gij$AA_GENO_PROB + Gij$AB_GENO_PROB +
             Gij$BA_GENO_PROB + Gij$BB_GENO_PROB),5) != 1, na.rm = T
    ) != 0
  ) {
    stop("AA_GENO_PROB + AB_GENO_PROB + BA_GENO_PROB + BB_GENO_PROB must equal 1 in all rows of Gij")
  }
  
  #Ensure that SAMPLE_ID column contains no 0s or NA
  if(
    sum(
      (is.na(Gij$SAMPLE_ID) |
       Gij$SAMPLE_ID == 0), na.rm = T
    ) != 0
  ) {
    stop("SAMPLE_ID in Gij cannot be NA or 0")
  } 
  
  Dij <- Gij
  colnames(Dij) <- c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB_G", "AB_GENO_PROB_G", "BA_GENO_PROB_G", "BB_GENO_PROB_G")
  
  #only retain data for genotype with maximum probability
  Dij$ROW_MAX <- apply(Dij[, c("AA_GENO_PROB_G", "AB_GENO_PROB_G", "BA_GENO_PROB_G", "BB_GENO_PROB_G")], 1, max) 
  Dij[!Dij[,"AA_GENO_PROB_G"] == Dij$ROW_MAX & !is.na(Dij[,"AA_GENO_PROB_G"]),"AA_GENO_PROB_G"] <- 0
  Dij[!Dij[,"AB_GENO_PROB_G"] == Dij$ROW_MAX & !is.na(Dij[,"AB_GENO_PROB_G"]),"AB_GENO_PROB_G"] <- 0
  Dij[!Dij[,"BA_GENO_PROB_G"] == Dij$ROW_MAX & !is.na(Dij[,"BA_GENO_PROB_G"]),"BA_GENO_PROB_G"] <- 0
  Dij[!Dij[,"BB_GENO_PROB_G"] == Dij$ROW_MAX & !is.na(Dij[,"BB_GENO_PROB_G"]),"BB_GENO_PROB_G"] <- 0
  
  #new columns (NA where GENO_PROB_G is NA)
  Dij[,"AA_GENO_PROB"] <- Dij[,"AA_GENO_PROB_G"] * 0
  Dij[,"AB_GENO_PROB"] <- Dij[,"AB_GENO_PROB_G"] * 0
  Dij[,"BA_GENO_PROB"] <- Dij[,"BA_GENO_PROB_G"] * 0
  Dij[,"BB_GENO_PROB"] <- Dij[,"BB_GENO_PROB_G"] * 0
  
  #Discrete genotypes from GENO_PROB_G if above threshold.indiv.  
  
  #If Gij is symetrical.
  if(identical(Dij[, "AB_GENO_PROB_G"],Dij[, "BA_GENO_PROB_G"])) {
    Dij[Dij[,"AA_GENO_PROB_G"] > threshold.indiv & !is.na(Dij[,"AA_GENO_PROB_G"]), "AA_GENO_PROB"] <- 1
    Dij[(Dij[,"AB_GENO_PROB_G"] + Dij[,"AB_GENO_PROB_G"]) > threshold.indiv & 
          !is.na(Dij[,"AB_GENO_PROB_G"]) & 
          !is.na(Dij[,"BA_GENO_PROB_G"]), "AB_GENO_PROB"] <- 1/2
    Dij[, "BA_GENO_PROB"] <- Dij[, "AB_GENO_PROB"] 
    Dij[Dij[,"BB_GENO_PROB_G"] > threshold.indiv & !is.na(Dij[,"BB_GENO_PROB_G"]), "BB_GENO_PROB"] <- 1
  }
  
  #If Gij is not symetrical.
  if(!identical(Dij[, "AB_GENO_PROB_G"],Dij[, "BA_GENO_PROB_G"])) {
    Dij[Dij[,"AA_GENO_PROB_G"] > threshold.indiv & !is.na(Dij[,"AA_GENO_PROB_G"]), "AA_GENO_PROB"] <- 1
    Dij[Dij[,"AB_GENO_PROB_G"] > threshold.indiv & !is.na(Dij[,"AB_GENO_PROB_G"]), "AB_GENO_PROB"] <- 1
    Dij[Dij[,"BA_GENO_PROB_G"] > threshold.indiv & !is.na(Dij[,"BA_GENO_PROB_G"]), "BA_GENO_PROB"] <- 1
    Dij[Dij[,"BB_GENO_PROB_G"] > threshold.indiv & !is.na(Dij[,"BB_GENO_PROB_G"]), "BB_GENO_PROB"] <- 1
  }
  
  Dij <- Dij[,c("SAMPLE_ID", "SNP_ID", "AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")]
  
  #probabilities must sum to 1
  row.sum <- rowSums(Dij[,c("AA_GENO_PROB", "AB_GENO_PROB", "BA_GENO_PROB", "BB_GENO_PROB")])
  Dij[row.sum != 1 & !is.na(row.sum),"AA_GENO_PROB"] <- NA
  Dij[row.sum != 1 & !is.na(row.sum),"AB_GENO_PROB"] <- NA
  Dij[row.sum != 1 & !is.na(row.sum),"BA_GENO_PROB"] <- NA
  Dij[row.sum != 1 & !is.na(row.sum),"BB_GENO_PROB"] <- NA
  
  #Equation 2 of Henshall et al. 2014 - T vector
  Dij$A_TRANS_PROB <- (Dij$AA_GENO_PROB * 2 + Dij$AB_GENO_PROB + Dij$BA_GENO_PROB)/2  #(sum row 1 + sum col 1) / 2
  Dij$B_TRANS_PROB <- (Dij$BB_GENO_PROB * 2 + Dij$AB_GENO_PROB + Dij$BA_GENO_PROB)/2  #(sum row 2 + sum col 2) / 2
  
  return(Dij)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

dkj.from.gkj.fun <- function(lambda.kj, threshold.pools) {
  
  #Description: Returns the Dij matrix for each individual and SNP according to the description on page 5 of 
  #Henshall et al. 2014 from the Gij matrices.  Discrete genotypes assigned (1 for homozygote, 0.5 for heterozygote)
  #to an element of Dij if the corresponding element of Gij greater than a threshold.pools.
  
  #Args##########################################
  
  #Gij: Data frame.  Output of Gij.fun 
  #          1.   SAMPLE_ID is the individual identifier (assumed to be offspring if not 
  #               listed as a sire or dam in 'fams' input)
  #          2.   SNP_ID is the SNP identifier
  #          3- Genotype probabilities (see the 
  #               top left of page 4 of Henshall et al. 2014
  
  #threshold.pools: Number.  Genotype called if AA_GENO_PROB, AB_GENO_PROB + BA_GENO_PROB or BB_GENO_PROB greater than this value
  
  #Returns##########################################
  # D.matrices:     Data frame 
  #              1.   SAMPLE_ID is the individual identifier 
  #              2.   SNP_ID is the SNP identifier
  #              3. Genotype probabilities: elements of the dkj matrix 
  
  print("Running dkj.from.gkj.fun")
  
  
  #  tmp <- gkj
  
  #get rho.inv
  #for(geno in rho.inv[,"GENOTYPE"]) {
  #   tmp[,geno] <- tmp[,geno] / rho.inv[rho.inv[,"GENOTYPE"]==geno,"RHO_INV"]
  #  }
  
  tmp <- lambda.kj[,-ncol(lambda.kj)]
  for(i in 3:ncol(tmp)) {
    colnames(tmp)[i] <- gsub("_LAMBDA","",colnames(tmp)[i])
  }
  tmp[,3:ncol(tmp)] <- tmp[,3:ncol(tmp)]/rowSums(tmp[,3:ncol(tmp)])
  
  n.in.pools <- nchar(colnames(tmp)[ncol(tmp)])/2
  
  rho.inv <- rho.inv.fun(n.in.pools)
  
  #get row maximums
  row.max <-  do.call(pmax, tmp[,3:ncol(tmp)])
  
  dkj <- tmp
  
  #Assign genotype if probability is the maximum in the row and greater than the threshold.pools
  for(geno in rho.inv[,"GENOTYPE"]) {
    dkj[tmp[,geno] == row.max &
          !is.na(tmp[,geno]) &
          row.max >= threshold.pools, geno] <- 1
    dkj[tmp[,geno] != row.max &
          !is.na(tmp[,geno]) &
          row.max >= threshold.pools, geno] <- 0
    dkj[row.max < threshold.pools |
          is.na(row.max), geno] <- NA
  }
  
  #probabilities must sum to 1
  row.sum <- round(rowSums(dkj[,!colnames(dkj) %in% c("SAMPLE_ID", "SNP_ID")]),10)
  dkj[row.sum != 1 & !is.na(row.sum),!colnames(dkj) %in% c("SAMPLE_ID", "SNP_ID")] <- NA
  
  #multiply by rho.inv because genotypes are unordered in dkj
  for(geno in rho.inv[,"GENOTYPE"]) {
    dkj[,geno] <- dkj[,geno] * rho.inv[rho.inv[,"GENOTYPE"]==geno,"RHO_INV"]
  }
  
  return(dkj)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

mismatches.fun <- function(dkj,
                           tclj,
                           miss.parent.count,
                           parent.combns,
                           fam.set.combns.by.pool) {
  
  #Exclusion method
  tmp <- duos.mismatches.fun(dkj = dkj,
                             tclj = tclj,
                             miss.parent.count = miss.parent.count,
                             fam.set.combns.by.pool =fam.set.combns.by.pool)
  duos.mismatches          <- tmp$duos.mismatches
  duos.mismatches.by.snp   <- tmp$duos.mismatches.by.snp
  mismatches.by.snp.sample <- tmp$mismatches.by.snp.sample
  mismatch.snp.count       <- tmp$mismatch.snp.count
  rm(tmp)
  
  most.like.parents <- most.like.parents.duos.mismatches.fun(duos.mismatches = duos.mismatches,
                                                             parent.combns = parent.combns)
  
  #Generate list where only common parents and families to all equally most likely parental combinations are retained
  
  most.like.parents.non.dup <- unique(most.like.parents[,-grep("ALT_",colnames(most.like.parents))]) #remove ALT parent combinations
  most.like.parents.non.dup <- unique(most.like.parents.non.dup[, !colnames(most.like.parents.non.dup) %in% c("PARENT_COMBN_ID", "FAM_COMBN_ID")])
  
  dup.samples <- unique(most.like.parents.non.dup[duplicated(most.like.parents.non.dup[,
                                                                                       c("SAMPLE_ID", "MISMATCHES", "SNP_COUNT", "MISMATCH_PROP", "MISMATCH_PROP_SE", "MISMATCH_PROP_Z")]),"SAMPLE_ID"])
  dup.most.like.parents <- most.like.parents.non.dup[most.like.parents.non.dup[,"SAMPLE_ID"] %in% dup.samples,]
  non.dup.most.like.parents <- most.like.parents.non.dup[!most.like.parents.non.dup[,"SAMPLE_ID"] %in% dup.samples,]
  tmp.dup.most.like.parents <- NULL
  
  for(samp in dup.samples) {
    tmp <- dup.most.like.parents[dup.most.like.parents[,"SAMPLE_ID"] == samp,]
    tmp.n <- nrow(tmp)
    
    tmp.parents <- tmp[,grep("PARENT_",colnames(most.like.parents.non.dup))]
    tmp.colnames.parents <- colnames(tmp.parents)
    tmp.parents <- matrix(unlist(tmp.parents), ncol = 1)
    tmp.parents[is.na(tmp.parents)] <- 0 #replace NA with 0
    tmp.parents <- as.data.frame(table(tmp.parents))
    tmp.parents[,"tmp.parents"] <- as.integer(levels(tmp.parents[,"tmp.parents"] ))
    tmp.parents[,"Freq"] <- tmp.parents[,"Freq"]/tmp.n
    tmp.parents[,"Freq"] <- floor(tmp.parents[,"Freq"])
    tmp.parents <- tmp.parents[tmp.parents[,"Freq"] > 0,]
    tmp.parents <- rep(tmp.parents[,1], tmp.parents[,2]) #list of parents common to all most likely combinations
    tmp.parents <- c(tmp.parents,rep(NA,length(tmp.colnames.parents) - length(tmp.parents) ))
    tmp.parents <- as.data.frame(t(tmp.parents))
    tmp.parents[tmp.parents == 0] <- NA #replace 0 with NA
    colnames(tmp.parents) <- tmp.colnames.parents
    
    tmp.fams <- tmp[,grep("FAMILY_ID_",colnames(most.like.parents.non.dup))]
    tmp.colnames.fams <- colnames(tmp.fams)
    tmp.fams <- matrix(unlist(tmp.fams), ncol = 1)
    tmp.fams[is.na(tmp.fams)] <- 0 #replace NA with 0
    tmp.fams <- as.data.frame(table(tmp.fams))
    tmp.fams[,"tmp.fams"] <- as.integer(levels(tmp.fams[,"tmp.fams"] ))
    tmp.fams[,"Freq"] <- tmp.fams[,"Freq"]/tmp.n
    tmp.fams[,"Freq"] <- floor(tmp.fams[,"Freq"])
    tmp.fams <- tmp.fams[tmp.fams[,"Freq"] > 0,]
    tmp.fams <- rep(tmp.fams[,1], tmp.fams[,2]) #list of families common to all most likely combinations
    tmp.fams <- c(tmp.fams,rep(NA,length(tmp.colnames.fams) - length(tmp.fams) ))
    tmp.fams <- as.data.frame(t(tmp.fams))
    tmp.fams[tmp.fams == 0] <- NA #replace 0 with NA
    colnames(tmp.fams) <- tmp.colnames.fams
    if(length(tmp.fams) == 0) {
      tmp.fams <- t(as.data.frame(rep(NA, length(tmp.parents)/2)))
      colnames(tmp.fams) <- paste("FAMILY_ID_", 1:(length(tmp.parents)/2), sep = "")
    }
    
    tmp <-  cbind(tmp[1,c("SAMPLE_ID", "MISMATCHES", "SNP_COUNT", "MISMATCH_PROP", "MISMATCH_PROP_SE", "MISMATCH_PROP_Z")],
                  tmp.fams, 
                  tmp.parents)
    tmp.dup.most.like.parents <- rbind(tmp.dup.most.like.parents,tmp)
    
    rm(tmp.colnames.parents, tmp.colnames.fams, tmp, tmp.n, tmp.fams, tmp.parents)
    
  }
  
  most.like.parents.non.dup <- rbind(tmp.dup.most.like.parents, non.dup.most.like.parents)
  rm(dup.most.like.parents, non.dup.most.like.parents, tmp.dup.most.like.parents)
  most.like.parents.non.dup <- most.like.parents.non.dup[order(most.like.parents.non.dup$SAMPLE_ID),]
  
  return(list(duos.mismatches = duos.mismatches,
              duos.mismatches.by.snp = duos.mismatches.by.snp,
              mismatches.by.snp.sample = mismatches.by.snp.sample,
              mismatch.snp.count = mismatch.snp.count,
              most.like.parents = most.like.parents,
              most.like.parents.non.dup = most.like.parents.non.dup))
}


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

duos.mismatches.fun <- function(dkj,
                                tclj,
                                miss.parent.count,
                                fam.set.combns.by.pool) {  
  
  #Computes count of mismatches by duo using the exclusion approach 
  
  print("Running duos.mismatches.fun")
  
  #identify missing SAMPLE_ID and SNP_ID combinations
  dkj$MISS_POOL <- rowSums(is.na(dkj[,3:ncol(dkj)])) != 0
  #  tclj <- merge(tclj, 
  #               miss.parent.count, 
  #               by = c("PARENT_COMBN_ID", "SNP_ID"), 
  #               all.x = TRUE)
  
  miss.parent.count$SNP_ID    <- as.character(miss.parent.count$SNP_ID)
  tclj$SNP_ID    <- as.character(tclj$SNP_ID)  
  
  #seed data frame
  duos.mismatches <- NULL
  mismatches.by.snp.sample <- NULL
  
  #Loop through pools 
  #Could merge Tfj with Dij but likely to run into memory issues
  
  pools <- unique(dkj$SAMPLE_ID)
  
  for(i in 1:length(pools)) {
    
    pool <- pools[i] #pools identifier    
    fam.set.combn <-fam.set.combns.by.pool[fam.set.combns.by.pool[,"SAMPLE_ID"] == pool,"FAM_SET_COMBN_ID"]
    tcj <- tclj[tclj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    
    tcj <- left_join(tcj, miss.parent.count, by = c("PARENT_COMBN_ID", "SNP_ID"))
    tcj$MISS_PARENT <- tcj$MISS_PARENT_COUNT > 0
    
    
    print(paste("Computing duo mismatches for pool", pool, "-", i, "of", length(pools)))
    
    mismatches <- merge(tcj,
                        dkj[dkj[,"SAMPLE_ID"] == pool,],
                        by = "SNP_ID",
                        all = TRUE,
                        suffixes = c("_TRANS","_POOLS")
    )
    
    #snp.count
    snp.count <- aggregate(!mismatches$MISS_PARENT & !mismatches$MISS_POOL,
                           by = list(mismatches$SAMPLE_ID,
                                     mismatches$PARENT_COMBN_ID), 
                           na.rm=T, 
                           FUN = "sum")   
    colnames(snp.count) <- c("SAMPLE_ID", "PARENT_COMBN_ID", "SNP_COUNT")
    
    #Equation 3 of Henshall et al. 2014
    mismatches[,"MISMATCHES"] <- rowSums(mismatches[,grepl("TRANS", colnames(mismatches))] * 
                                           mismatches[,grepl("POOLS", colnames(mismatches))])
    mismatches[,"MISMATCHES"] <- as.integer(mismatches[,"MISMATCHES"] == 0)
    
    mismatches[mismatches[,"MISS_PARENT"] == TRUE, "MISMATCHES"] <- NA
    mismatches[mismatches[,"MISS_POOL"]   == TRUE, "MISMATCHES"] <- NA
    
    snp.mismatches <- aggregate(mismatches[,"MISMATCHES"], by = list(mismatches$SAMPLE_ID,
                                                                     mismatches$PARENT_COMBN_ID), 
                                na.rm=T, FUN = "sum")     
    colnames(snp.mismatches) <- c("SAMPLE_ID", "PARENT_COMBN_ID", "MISMATCHES")
    
    #   snp.mismatches <- merge(snp.mismatches,snp.count, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.x = TRUE)
    snp.mismatches <- left_join(snp.mismatches, snp.count, by = c("SAMPLE_ID", "PARENT_COMBN_ID"))
    
    #Calculate proportion of SNPs that are mismatched
    snp.mismatches$MISMATCH_PROP <- snp.mismatches$MISMATCHES / snp.mismatches$SNP_COUNT
    snp.mismatches$MISMATCH_PROP_SE <- sqrt((snp.mismatches$MISMATCH_PROP*(1-snp.mismatches$MISMATCH_PROP)) / snp.mismatches$SNP_COUNT)
    
    #identify most likely to identify bad SNP
    
    #Obtain most likely combn (i.e. min MISMATCH_PROP) for each pool
    
    #if no SNP for which genotypes were assigned in pool
    if(max(snp.mismatches$SNP_COUNT) == 0) {
      mismatches.by.snp <- mismatches
      mismatches.by.snp <- mismatches.by.snp[mismatches.by.snp[,"SAMPLE_ID"] == pool,]
      mismatches.by.snp[,!colnames(mismatches.by.snp) %in% c("SNP_ID", "SAMPLE_ID")] <- NA
      mismatches.by.snp <- unique(mismatches.by.snp)
      
      duos.mismatches <- rbind(duos.mismatches,snp.mismatches)
      mismatches.by.snp.sample <- rbind(mismatches.by.snp.sample,mismatches.by.snp)
      print(paste("WARNING: Exclusion approach to pedigree assignment is not functioning appropriately as there are no SNP for which genotypes were assigned in pool", pool, ". Consider using changing discrete.method to equal \'assigned.genos\' or reducing the value of threshold.pools"))
    } else {
      
      most.like <-  min(snp.mismatches$MISMATCH_PROP, na.rm = TRUE) 
      most.like <- snp.mismatches[snp.mismatches[,"MISMATCH_PROP"] == most.like,]
      
      #if there are multiple combns then go with the one with the most snp
      if(nrow(most.like) > 1) {
        tmp.most.like <-  aggregate(most.like$SNP_COUNT, 
                                    by = list(most.like$SAMPLE_ID), na.rm=T, FUN = "max") 
        colnames(tmp.most.like) <- c("SAMPLE_ID", "SNP_COUNT")
        
        #   most.like <- merge(tmp.most.like, most.like, by = c("SAMPLE_ID", "SNP_COUNT"), all.x = TRUE)
        most.like <- left_join(tmp.most.like, most.like, by = c("SAMPLE_ID", "SNP_COUNT"))
        
        most.like <- most.like[,c("SAMPLE_ID",	"PARENT_COMBN_ID", "MISMATCHES",	"SNP_COUNT",	
                                  "MISMATCH_PROP",	"MISMATCH_PROP_SE")]
        #remove duplicates at random
        most.like <- most.like[order(runif(nrow(most.like))),]
        most.like <- most.like[most.like[,"PARENT_COMBN_ID"] == most.like[1,"PARENT_COMBN_ID"],]
      }
      
      #add mismatched SNP in most.like to total counts of bad snp
      
      mismatches.by.snp <- mismatches[mismatches[,"PARENT_COMBN_ID"] == most.like[,"PARENT_COMBN_ID"],]
    }
    duos.mismatches <- rbind(duos.mismatches,snp.mismatches)
    mismatches.by.snp.sample <- rbind(mismatches.by.snp.sample,mismatches.by.snp)
    
  } #END for(i in 1:length(pools)) {
  
  duos.mismatches$MISMATCH_PROP_Z <- duos.mismatches$MISMATCH_PROP / duos.mismatches$MISMATCH_PROP_SE 
  duos.mismatches[duos.mismatches[,"MISMATCH_PROP"] == 0 &
                    !is.na(duos.mismatches[,"MISMATCH_PROP"]),  
                  "MISMATCH_PROP_Z"] <- 0
  
  duos.mismatches <- duos.mismatches[order(duos.mismatches$MISMATCH_PROP),]
  duos.mismatches <- duos.mismatches[order(duos.mismatches$SAMPLE_ID),]
  
  #Count mismatches
  
  mismatch.snp.count <- aggregate(mismatches.by.snp.sample[,"MISMATCHES"] == 1, 
                                  by = list(mismatches.by.snp.sample$SNP_ID), 
                                  na.rm=T, FUN = "sum")   
  colnames(mismatch.snp.count) <- c("SNP_ID", "MISMATCH_COUNT")
  
  match.snp.count <- aggregate(mismatches.by.snp.sample[,"MISMATCHES"] == 0, 
                               by = list(mismatches.by.snp.sample$SNP_ID), 
                               na.rm=T, FUN = "sum")   
  colnames(match.snp.count) <- c("SNP_ID", "MATCH_COUNT")
  
  miss.snp.count <- aggregate(is.na(mismatches.by.snp.sample[,"MISMATCHES"]), 
                              by = list(mismatches.by.snp.sample$SNP_ID), 
                              na.rm=T, FUN = "sum")   
  colnames(miss.snp.count) <- c("SNP_ID", "MISSING_COUNT")
  
  #  mismatch.snp.count <- merge(mismatch.snp.count, match.snp.count, by = "SNP_ID", all = TRUE)  
  mismatch.snp.count$SNP_ID <- as.character(mismatch.snp.count$SNP_ID)
  match.snp.count$SNP_ID    <- as.character(match.snp.count$SNP_ID)
  mismatch.snp.count        <- inner_join(mismatch.snp.count, match.snp.count, by = "SNP_ID")
  
  #mismatch.snp.count <- merge(mismatch.snp.count, miss.snp.count, by = "SNP_ID", all = TRUE)  
  mismatch.snp.count$SNP_ID <- as.character(mismatch.snp.count$SNP_ID)
  miss.snp.count$SNP_ID    <- as.character(miss.snp.count$SNP_ID)
  mismatch.snp.count        <- inner_join(mismatch.snp.count, miss.snp.count, by = "SNP_ID")
  
  mismatch.snp.count <- mismatch.snp.count[order(mismatch.snp.count$MISSING_COUNT, decreasing = TRUE),]
  mismatch.snp.count <- mismatch.snp.count[order(mismatch.snp.count$MISMATCH_COUNT, decreasing = TRUE),]
  
  return(list(duos.mismatches = duos.mismatches, 
              duos.mismatches.by.snp = mismatches, 
              mismatches.by.snp.sample = mismatches.by.snp.sample,
              mismatch.snp.count = mismatch.snp.count))
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export
#' 
most.like.parents.duos.mismatches.fun <- function(duos.mismatches,
                                                  parent.combns) {
  
  #Description: Identifies the family that each inidividual is most likely to belong to.
  
  #  Args:
  
  #  Returns: data frame
  #          1.  SAMPLE_ID 
  #          2.  PARENT_COMBN_ID 
  #              FAM_COMBN_ID
  #          3.  MISS_PARENT_SNP_DATA_PROP 
  #          4.  MISS_POOL_SNP_DATA_PROP 
  #          5.  NO_MISS_PARENT_OR_POOL_PROP 
  #              PARENT_1 PARENT_2 PARENT_3 PARENT_4 ...
  #          6.  SNP_COUNT: Count of SNP  (only includes SNP for which genotypes for the sire, dam and pool were present)
  #          7.  MISMATCHES: Mismatched snp (only includes SNP for which genotypes for the sire, dam and pool were present)
  #          8.  MISMATCH_PROP: MISMATCHES / SNP_COUNT
  #          9.  MISMATCH_PROP_SE: Standard error of MISMATCH_PROP (where SE = sqrt [ p(1 - p) / n ] )
  #          10. MISMATCH_PROP_Z
  #              ALT_PARENT_COMBN_ID 
  #              ALT_FAM_COMBN_ID 
  #              ALT_PARENT_1 ALT_PARENT_2 ALT_PARENT_3 ALT_PARENT_4 ...
  #          11. ALT_SNP_COUNT: Count of SNP  (only includes SNP for which genotypes for the sire, dam and pool were present)
  #          12. ALT_MISMATCHES: Mismatched snp (only includes SNP for which genotypes for the sire, dam and pool were present)
  #          13. ALT_MISMATCH_PROP: MISMATCHES / SNP_COUNT
  #          14. ALT_MISMATCH_PROP_SE: Standard error of MISMATCH_PROP (where SE = sqrt [ p(1 - p) / n ] )
  #          15. ALT_MISMATCH_PROP_Z  
  
  print("Running most.like.parents.duos.mismatches.fun")
  
  #remove FAM_COMBN_ID containing duplicated FAM_COMBN_ID
  tmp <- parent.combns[duplicated(parent.combns[,"PARENT_COMBN_ID"]),
                       "PARENT_COMBN_ID"]
  
  if(length(tmp) != 0) {
    parent.combns[parent.combns[,"PARENT_COMBN_ID"] %in% tmp,grep("FAM",colnames(parent.combns))] <- NA
  }
  rm(tmp)
  parent.combns.unambiguous <- parent.combns
  
  #Obtain most likely combn (i.e. min MISMATCH_PROP) for each pool
  most.like           <-  aggregate(duos.mismatches$MISMATCH_PROP, 
                                    by = list(duos.mismatches$SAMPLE_ID), na.rm=T, FUN = "min")   
  colnames(most.like) <- c("SAMPLE_ID", "MISMATCH_PROP")
  
  #  most.like <-  merge(most.like, duos.mismatches, by = c("SAMPLE_ID", "MISMATCH_PROP"), all.x = TRUE)
  most.like <- left_join(most.like, duos.mismatches, by = c("SAMPLE_ID", "MISMATCH_PROP"))
  
  #if there are multiple families then go with the one with the most snp
  tmp.most.like <-  aggregate(most.like$SNP_COUNT, 
                              by = list(most.like$SAMPLE_ID), na.rm=T, FUN = "max") 
  colnames(tmp.most.like) <- c("SAMPLE_ID", "SNP_COUNT")
  
  # most.like <-  merge(tmp.most.like, most.like, by = c("SAMPLE_ID", "SNP_COUNT"), all.x = TRUE)
  most.like               <- left_join(tmp.most.like, most.like, by = c("SAMPLE_ID", "SNP_COUNT"))
  
  most.like <- most.like[,c("SAMPLE_ID",	"PARENT_COMBN_ID", "MISMATCHES",	"SNP_COUNT",	
                            "MISMATCH_PROP",	"MISMATCH_PROP_SE",	"MISMATCH_PROP_Z")]
  # most.like <-  merge(most.like, parent.combns.unambiguous, by = c("PARENT_COMBN_ID"), all.x = TRUE)
  most.like <- left_join(most.like, parent.combns.unambiguous, by = "PARENT_COMBN_ID")
  
  most.like <- unique(most.like)
  
  #Remove most likely trio to obtain second most likely
  remove <- most.like[,c("SAMPLE_ID", "PARENT_COMBN_ID")]
  remove$REMOVE <- TRUE
  # duos.mismatches.removed <- merge(duos.mismatches, remove, by = c("SAMPLE_ID", "PARENT_COMBN_ID"), all.x = TRUE)
  duos.mismatches.removed   <- left_join(duos.mismatches, remove, by = c("SAMPLE_ID", "PARENT_COMBN_ID"))
  duos.mismatches.removed   <- duos.mismatches.removed[is.na(duos.mismatches.removed$REMOVE),]
  
  #Obtain second most likely family (i.e. maximum LOD) for each pool
  second.most.like <-  aggregate(duos.mismatches.removed$MISMATCH_PROP, 
                                 by = list(duos.mismatches.removed$SAMPLE_ID), na.rm=T, FUN = "min")   
  colnames(second.most.like) <- c("SAMPLE_ID", "MISMATCH_PROP")
  
  #  second.most.like <-  merge(second.most.like, duos.mismatches.removed, 
  #                             by = c("SAMPLE_ID", "MISMATCH_PROP"), all.x = TRUE)
  second.most.like   <- left_join(second.most.like, duos.mismatches.removed, by = c("SAMPLE_ID", "MISMATCH_PROP"))
  
  tmp.second.most.like <-  aggregate(second.most.like$SNP_COUNT, 
                                     by = list(second.most.like$SAMPLE_ID), na.rm=T, FUN = "max") 
  colnames(tmp.second.most.like) <- c("SAMPLE_ID", "SNP_COUNT")
  
  # second.most.like <-  merge(tmp.second.most.like, second.most.like, by = c("SAMPLE_ID", "SNP_COUNT"), all.x = TRUE)
  second.most.like <- left_join(tmp.second.most.like, second.most.like, by = c("SAMPLE_ID", "SNP_COUNT"))
  
  #  second.most.like <- merge(second.most.like, 
  #                             parent.combns.unambiguous[,c("PARENT_COMBN_ID", "FAM_COMBN_ID")], 
  #                             by = c("PARENT_COMBN_ID"), all.x = TRUE)
  second.most.like <- left_join(second.most.like, 
                                parent.combns.unambiguous[,c("PARENT_COMBN_ID", "FAM_COMBN_ID")], 
                                by = "PARENT_COMBN_ID")
  
  second.most.like <- second.most.like[,c("SAMPLE_ID",	"PARENT_COMBN_ID",	"FAM_COMBN_ID", "MISMATCHES",	"SNP_COUNT",	
                                          "MISMATCH_PROP",	"MISMATCH_PROP_SE",	"MISMATCH_PROP_Z")]
  
  colnames(second.most.like) <- c("SAMPLE_ID",	"ALT_PARENT_COMBN_ID",	"ALT_FAM_COMBN_ID","ALT_MISMATCHES",	"ALT_SNP_COUNT",	
                                  "ALT_MISMATCH_PROP",	"ALT_MISMATCH_PROP_SE",	"ALT_MISMATCH_PROP_Z")
  second.most.like <- unique(second.most.like)
  
  #Merge most likely and second most likely 
  #  most.like <- merge(most.like, second.most.like, by = "SAMPLE_ID", all = TRUE)
  most.like <- inner_join(most.like, second.most.like, by = "SAMPLE_ID")
  
  return(most.like)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

ls.fun <- function(fams,
                   fam.set.combns,
                   fam.set.combns.by.pool,
                   Gij,
                   flj,
                   snp.dat.pools,
                   snp.param.indiv,
                   min.intensity,
                   beta.min.ss){
  
  #retain only required columns
  snp.dat.pools <- snp.dat.pools[,c("SNP_ID", "SAMPLE_ID", "INTENSITY_A", "INTENSITY_B")]
  
  #get parent.combns to generate tclj
  tmp.fam.set.combns <- fam.set.combns
  tmp.fam.set.combns$FAM_SET_ID=1 #Assume n.in.pools = 1 and there is only one family set to get appropriate 'parent.combns' (i.e. 1 row per family)
  tmp <- parent.combns.fun(fams = fams,
                           n.in.pools = 1, 
                           fam.set.combns = tmp.fam.set.combns)
  parent.combns.by.fam.set.combn <- tmp$parent.combns.by.fam.set.combn
  parent.combns <- tmp$parent.combns 
  rm(tmp.fam.set.combns)
  
  tclj <- NULL
  for(fam.set.combn in unique(flj$FAM_SET_COMBN_ID)) {
    
    tmp.fj <- flj[flj[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    tmp.parent.combns <- parent.combns.by.fam.set.combn[parent.combns.by.fam.set.combn[,"FAM_SET_COMBN_ID"] == fam.set.combn,] 
    
    tmp.parents <- tmp.parent.combns[,!colnames(tmp.parent.combns) %in% c("FAM_SET_COMBN_ID", "PARENT_COMBN_ID")]
    tmp.parents <- unique(as.vector(as.matrix(tmp.parents)))
    
    tmp.tij <- Gij[Gij[,"SAMPLE_ID"] %in% tmp.parents, ]
    tmp.tclj <- tcj.fun(tij = tmp.tij, 
                        parent.combns = tmp.parent.combns, 
                        fj = tmp.fj)
    tmp.tclj$FAM_SET_COMBN_ID <- fam.set.combn
    tclj <- rbind(tclj, tmp.tclj)
  }
  rm(tmp.fj, tmp.parent.combns, tmp.parents, tmp.tij, tmp.tclj)
  
  #Get B_TRANS_PROB
  tclj$B_TRANS_PROB <- tclj$AB + tclj$BB #note that tclj$AB was divided by two in tcj.fun
  tclj <- tclj[,c("SNP_ID", "PARENT_COMBN_ID", "MISS_PARENT_COUNT", "FAM_SET_COMBN_ID", "B_TRANS_PROB")]
  
  fkj.and.weight <- fkj.and.weight.fun(snp.dat.pools = snp.dat.pools, 
                                       snp.param.indiv = snp.param.indiv, 
                                       min.intensity = min.intensity)
   
  Xl.mat <- X.mat.fun(tclj.ls = tclj, parent.combns = parent.combns)
  
  beta <- NULL
  #Loop through samples to get beta 
  for(samp in unique(fam.set.combns.by.pool[,"SAMPLE_ID"])) {
    
    fam.set.combn <-fam.set.combns.by.pool[fam.set.combns.by.pool[,"SAMPLE_ID"] == samp, "FAM_SET_COMBN_ID"]
    
    X.mat <- Xl.mat[[fam.set.combn]]
    #fkj.and.weight with samples from fam.set.combn only
    tmp.fkj.and.weight <- fkj.and.weight[fkj.and.weight[,"SAMPLE_ID"] == samp,]
    
    if(!beta.min.ss) {
      tmp.beta <- beta.fun(X.mat = X.mat,
                           fkj.and.weight = tmp.fkj.and.weight,
                           fams = fams,
                           miss.x.and.y = 0.5, 
                           miss.w = 0)
      tmp.beta$BETA_MIN_SS <- NA
    } else {
      tmp.beta <- beta.min.ss.fun(X.mat = X.mat,
                                  fkj.and.weight = tmp.fkj.and.weight,
                                  fam.set.combns = fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == fam.set.combn,],
                                  fams = fams)
    }
    
    #get BETA_HAT_CONSTRAINED
    tmp.fam.set.combns <- fam.set.combns.by.pool[fam.set.combns.by.pool[,"SAMPLE_ID"] == samp, "FAM_SET_COMBN_ID"]
    tmp.fam.set.combns <- fam.set.combns[fam.set.combns[,"FAM_SET_COMBN_ID"] == tmp.fam.set.combns, c("FAM_SET_ID", "FAMILY_ID")]
    beta.hat.constrained <- NULL
    n.in.pool <- length(unique(tmp.fam.set.combns[,"FAM_SET_ID"]))
    tmp <- tmp.beta
    for(fam.set in unique(tmp.fam.set.combns[,"FAM_SET_ID"])) {
      tmp.beta.hat.constrained <- tmp[tmp[,"FAMILY_ID"] %in% 
                                        tmp.fam.set.combns[tmp.fam.set.combns[,"FAM_SET_ID"] == fam.set,"FAMILY_ID"],]
      tmp.beta.hat.constrained <- tmp.beta.hat.constrained[tmp.beta.hat.constrained[,"BETA_HAT"] == 
                                                             max(tmp.beta.hat.constrained[,"BETA_HAT"]),"FAMILY_ID"] #identify familiy with maximum beta
      
      tmp[tmp[,"FAMILY_ID"] == tmp.beta.hat.constrained,"BETA_HAT"] <- 
        tmp[tmp[,"FAMILY_ID"] ==  tmp.beta.hat.constrained,"BETA_HAT"] - 1/n.in.pool
      beta.hat.constrained <- c(beta.hat.constrained,tmp.beta.hat.constrained)
    }
    beta.hat.constrained <- data.frame(FAMILY_ID = beta.hat.constrained,
                                       BETA_HAT_CONSTRAINED = 1/n.in.pool)
    beta.hat.constrained <- aggregate(BETA_HAT_CONSTRAINED ~ FAMILY_ID, data = beta.hat.constrained, FUN = sum) #sum BETA_HAT_CONSTRAINED by family
    tmp.beta <- merge(tmp.beta, beta.hat.constrained, by = "FAMILY_ID", all.x = TRUE)
    tmp.beta[is.na(tmp.beta[,"BETA_HAT_CONSTRAINED"]),"BETA_HAT_CONSTRAINED"] <- 0
    rm(beta.hat.constrained, tmp.beta.hat.constrained)
    
    beta <- rbind(beta,tmp.beta)
    rm(tmp.fkj.and.weight, tmp.beta)
  }
  
  beta <- beta[,c("SAMPLE_ID", "SIRE_ID", "DAM_ID", "FAMILY_ID", "BETA_STAR",  "BETA_HAT", "BETA_HAT_CONSTRAINED", "BETA_MIN_SS")]
  
  return(list(tclj.ls = tclj,
              fkj.and.weight = fkj.and.weight,
              Xl.mat = Xl.mat,
              beta = beta
  ))
  
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

fkj.and.weight.fun <- function(snp.dat.pools, 
                               snp.param.indiv, 
                               min.intensity = 0) {
  
  #Estimates B allele frequencies for pools.  According to the method
  #of Henshall et al. 2014.  Note that "Unlike
  #in the approach used to estimate B allele Frequencies in [24], where allelic
  #proprtions were scaled such that values p(SAMPLE_ID) - 00, 0.5 and 1.0 corrensponded 
  #to the means AA, AB and BB, respectively, we did not scale the p(SAMPLE_ID) values"
  #page 3 of Henshall et al. 2014
  
  #NOTE: Assumed ab.FREQ = BA.FREQ in genotype probability matrix (i.e. it is symetrical)
  
  #Required functions:
  # pij.fun
  
  #Args##########################################
  # snp.dat.pools: Data frame (pooled samples only) containing relevant fields from the Genotype.Intensity 
  #tab of the corresponding GenotypeIntensity.xls file outputted from Sequenom's Typer software 
  #(Sequenom 2006).  Equivalent outputs from other platforms could also be used.  
  #              1. SAMPLE_ID  is the pool sample identifier 
  #              2. SNP_ID  is the SNP identifier,
  #              3. INTENSITY_A   is the area/intensity for allele A, 
  #              4. INTENSITY_B   is the area/intensity for allele B
  
  # snp.param.indiv: Data frame (see "Estimation of SNP specific parameters, page 3 of Henshall et al. 2014).
  # Data frame containing relevant fields from the output of snp.gen.param.fun from individual animal data (e.g. parents).
  #              1. SNP_ID    is the SNP identifier, 
  #              2. MEAN_P_AA is the mean of allelic proportion (homozygous allele A), 
  #              3. MEAN_P_AB is the mean of allelic proportion (heterozygous), 
  #              4. MEAN_P_BB is the mean of allelic proportion (homozygous allele B), 
  #              5. WELCH_A   is the welch statistics for the intervals of mean of 
  #                           allelic proportion AA to AB (WELCH_A) 
  #              6. WELCH_B   is the welch statistics for the intervals of mean of 
  #                           allelic proportion AB to BB 
  
  # min.intensity      Number used in pij.fun. If sqrt((snp.dat.pools$INTENSITY_A)^2 +
  #              (snp.dat.pools$INTENSITY_B)^2) less than this value
  #              then set allelic proportion to missing (see end of page 3 of Henshall et al 2014).
  #              Essentially removes observations that fall into the lower left of INTENSITY_A
  #              by INTENSITY_B scatter plot.
  
  
  #Returns##########################################
  # fkj.and.weight:   Data frame
  #              1. SNP_ID    is the SNP identifier
  #              2. SAMPLE_ID   is the pool identifier 
  #              3. MEAN_P_AA is the mean of allelic proportion (homozygous allele A) from snp.param.indiv
  #              4. MEAN_P_AB is the mean of allelic proportion (heterozygous) from snp.param.indiv
  #              5. MEAN_P_BB is the mean of allelic proportion (homozygous allele B) from snp.param.indiv
  #              6. ALLELIC_PROP_POOL is the allelic proportion for pool (Equation 1 of Henshall et al. 2014)
  #              7. FREQ_POOL is the estimated frequency of the B allele of the SNP in the pool (y vector of page 6 of Henshall et al. 2014)
  #              8. FREQ_POOL_ERROR_WT  is the relevant Welch statistic squared (y vector of page 6 of Henshall et al. 2014)
  
  print("Running fkj.and.weight.fun")
  
  if(sum(c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B") %in% colnames(snp.dat.pools)) != 4) {
    stop("snp.dat.pools input must be a data frame containing the following headings: SAMPLE_ID, SNP_ID, INTENSITY_A, INTENSITY_B")
  }
  
  if(sum(c("SNP_ID", "MEAN_P_AA", "MEAN_P_AB", "MEAN_P_BB", "WELCH_A", "WELCH_B") %in% colnames(snp.param.indiv)) != 6) {
    stop("snp.param.indiv input must be a data frame containing the following headings: SNP_ID, MEAN_P_AA, MEAN_P_AB, MEAN_P_BB, WELCH_A, WELCH_B")
  }
  
  #Name columns
  snp.dat.pools$SAMPLE_ID   <- as.integer(snp.dat.pools$SAMPLE_ID)
  snp.dat.pools$SNP_ID    <- as.character(snp.dat.pools$SNP_ID)
  snp.dat.pools$INTENSITY_A    <- as.numeric(snp.dat.pools$INTENSITY_A)
  snp.dat.pools$INTENSITY_B    <- as.numeric(snp.dat.pools$INTENSITY_B)
  snp.dat.pools <- snp.dat.pools[,c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")]
  
  snp.param.indiv$SNP_ID      <- as.character(snp.param.indiv$SNP_ID)
  snp.param.indiv$MEAN_P_AA     <- as.numeric(snp.param.indiv$MEAN_P_AA)
  snp.param.indiv$MEAN_P_AB     <- as.numeric(snp.param.indiv$MEAN_P_AB)
  snp.param.indiv$MEAN_P_BB     <- as.numeric(snp.param.indiv$MEAN_P_BB)
  snp.param.indiv$WELCH_A     <- as.numeric(snp.param.indiv$WELCH_A)
  snp.param.indiv$WELCH_B     <- as.numeric(snp.param.indiv$WELCH_B)
  snp.param.indiv <- snp.param.indiv[,c("SNP_ID", "MEAN_P_AA", "MEAN_P_AB", "MEAN_P_BB", "WELCH_A", "WELCH_B")]
  
  # Check the list of SNPs the same in input files
  if(sum(!unique(snp.dat.pools[,"SNP_ID"]) %in% unique(snp.param.indiv[,"SNP_ID"]))>0 &
     sum(!unique(snp.param.indiv[,"SNP_ID"]) %in% unique(snp.dat.pools[,"SNP_ID"]))>0) {
    print("SNP identifiers do not match in snp.dat.pools and snp.param.indiv")
    stop()
  }
  
  #Get allelic proportion
  snp.allelic.prop <- pij.fun(snp.dat.indiv = snp.dat.pools, min.intensity = min.intensity)
  print("Still running fkj.and.weight.fun")
  
  #Rename columns
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "SAMPLE_ID"]               <- "SAMPLE_ID"
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "ALLELIC_PROP"]     <- "ALLELIC_PROP_POOL"  
  colnames(snp.allelic.prop)[colnames(snp.allelic.prop) == "INTENSITY"]        <- "INTENSITY_POOL"  
  
  #Merge data using cbind
  if(identical(snp.dat.pools$SAMPLE_ID,snp.allelic.prop$SAMPLE_ID) &
     identical(snp.dat.pools$SNP_ID,snp.allelic.prop$SNP_ID)) {
    snp.dat.pools <- cbind(snp.dat.pools, snp.allelic.prop[,!colnames(snp.allelic.prop) %in% c("SAMPLE_ID", "SNP_ID", "INTENSITY_A", "INTENSITY_B")])
  } else {
    stop("SAMPLE_ID and SNP_ID columns of pij.fun output do not match those of snp.dat.pools.  Not sure why.")
  }
  
  #Merge data
  # fkj.and.weight <- merge(snp.dat.pools, snp.param.indiv, by = "SNP_ID", all.x = TRUE)
  snp.dat.pools$SNP_ID    <- as.character(snp.dat.pools$SNP_ID)
  snp.param.indiv$SNP_ID    <- as.character(snp.param.indiv$SNP_ID)
  fkj.and.weight <- left_join(snp.dat.pools, snp.param.indiv, by = "SNP_ID")
  
  #Compute allele frequencies (see rhs of page 5 of Henshall et al. 2014)
  fkj.and.weight$FREQ_POOL <- NA
  
  #ALLELIC_PROP_POOL <= MEAN_P_AA
  is.true <- fkj.and.weight[,"ALLELIC_PROP_POOL"] <= fkj.and.weight[,"MEAN_P_AA"]
  is.true[is.na(is.true)] <- FALSE
  
  fkj.and.weight[is.true,"FREQ_POOL"] <- 0
  fkj.and.weight[is.true,"FREQ_ERROR"]      <- fkj.and.weight[is.true,"WELCH_A"]^2 #WELCH_A squared
  
  #MEAN_P_AA < ALLELIC_PROP_POOL < MEAN_P_AB
  is.true <- fkj.and.weight[,"MEAN_P_AA"] < fkj.and.weight[,"ALLELIC_PROP_POOL"] &
    fkj.and.weight[,"ALLELIC_PROP_POOL"] < fkj.and.weight[,"MEAN_P_AB"] 
  is.true[is.na(is.true)] <- FALSE
  
  fkj.and.weight[is.true,"FREQ_POOL"] <- 0.5 *
    ((fkj.and.weight[is.true,"ALLELIC_PROP_POOL"] - fkj.and.weight[is.true,"MEAN_P_AA"]) / 
       (fkj.and.weight[is.true,"MEAN_P_AB"] - fkj.and.weight[is.true,"MEAN_P_AA"]))
  fkj.and.weight[is.true,"FREQ_ERROR"]      <- fkj.and.weight[is.true,"WELCH_A"]^2 #WELCH_A squared
  
  #MEAN_P_AB < ALLELIC_PROP_POOL < MEAN_P_BB
  is.true <- fkj.and.weight[,"MEAN_P_AB"] < fkj.and.weight[,"ALLELIC_PROP_POOL"] &
    fkj.and.weight[,"ALLELIC_PROP_POOL"] < fkj.and.weight[,"MEAN_P_BB"]
  is.true[is.na(is.true)] <- FALSE
  
  fkj.and.weight[is.true,"FREQ_POOL"] <- 0.5 + 0.5 *
    ((fkj.and.weight[is.true,"ALLELIC_PROP_POOL"] - fkj.and.weight[is.true,"MEAN_P_AB"]) / 
       (fkj.and.weight[is.true,"MEAN_P_BB"] - fkj.and.weight[is.true,"MEAN_P_AB"])) 
  fkj.and.weight[is.true,"FREQ_ERROR"]      <- fkj.and.weight[is.true,"WELCH_B"]^2 #WELCH_B squared
  
  #ALLELIC_PROP_POOL >= MEAN_P_BB
  is.true <- fkj.and.weight[,"ALLELIC_PROP_POOL"] >= fkj.and.weight[,"MEAN_P_BB"]
  is.true[is.na(is.true)] <- FALSE
  
  fkj.and.weight[is.true,"FREQ_POOL"] <- 1
  fkj.and.weight[is.true,"FREQ_ERROR"]      <- fkj.and.weight[is.true,"WELCH_B"]^2 #WELCH_B squared
  rm(is.true)
  
  #If ALLELIC_PROP_POOL, MEAN_P_AA, MEAN_P_AB or MEAN_P_BB = NA, then make FREQ_POOL = NA
  fkj.and.weight[is.na(fkj.and.weight[,"ALLELIC_PROP_POOL"]) |
                   is.na(fkj.and.weight[,"MEAN_P_AA"]) |
                   is.na(fkj.and.weight[,"MEAN_P_AB"]) |
                   is.na(fkj.and.weight[,"MEAN_P_BB"]),"FREQ_POOL"] <- NA
  
  #Retain relevant columns only
  fkj.and.weight <- fkj.and.weight[,c("SNP_ID", "SAMPLE_ID", "MEAN_P_AA", "MEAN_P_AB", "MEAN_P_BB",
                                      "ALLELIC_PROP_POOL", "FREQ_POOL", "FREQ_ERROR", "SAMPLE_ID")]
  colnames(fkj.and.weight)[colnames(fkj.and.weight) == "FREQ_ERROR"] <- "FREQ_POOL_ERROR_WT" 
  
  #Order by pool identifier
  fkj.and.weight <- fkj.and.weight[order(fkj.and.weight$SAMPLE_ID, decreasing = FALSE), ]   #help("order")
  
  #if FREQ_POOL is NA then make FREQ_POOL_ERROR_WT NA
  fkj.and.weight[is.na(fkj.and.weight[,"FREQ_POOL"]),"FREQ_POOL_ERROR_WT"] <- NA
  
  fkj.and.weight <- fkj.and.weight[,colnames(fkj.and.weight) != "SAMPLE_ID.1"]
  
  return(fkj.and.weight)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

X.mat.fun <- function(tclj.ls, parent.combns, miss.x.and.y = 0.5) {
  #Generate the X matrix, as described in Henshall et al. (2014)   
  # miss.x.and.y Number. Value to replace missing values in the X matrix and y vector (Henshall et al. 2014 page 6).  
  
  print("Running X.mat.fun")
  
 # if("reshape2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("reshape2")} 
  library(reshape2)
  
  if(sum(c("SNP_ID", "PARENT_COMBN_ID", "MISS_PARENT_COUNT", "B_TRANS_PROB") %in% colnames(tclj.ls)) != 4) {
    stop("tclj.ls input must be a data frame containing the following headings: SNP_ID, PARENT_COMBN_ID, MISS_PARENT_COUNT, B_TRANS_PROB")
  }
  
  fam.set.combns <- unique(tclj.ls[,"FAM_SET_COMBN_ID"])
  
  X.mat <- NULL
  
  for(fam.set.combn in fam.set.combns) {
    
    tcj <- tclj.ls[tclj.ls[,"FAM_SET_COMBN_ID"] == fam.set.combn,]
    tcj <- tcj[,colnames(tcj) != "FAM_SET_COMBN_ID"]
    
    #  tcj <- merge(tcj, parent.combns, by = "PARENT_COMBN_ID", all.x = TRUE)
    tcj <- left_join(tcj, parent.combns, by = "PARENT_COMBN_ID")
    
    colnames(tcj)[colnames(tcj) == "FAMILY_ID_1"] <- "FAMILY_ID"
    
    #Get B_TRANS_PROB
    # tcj$B_TRANS_PROB <- tcj$AB/2  + tcj$BB
    
    #Replace missing values
    #If missing they must be missing in fjfj and snp.param.indiv 
    tcj[tcj$MISS_PARENT_COUNT == 2,"B_TRANS_PROB"] <- miss.x.and.y #if both sire and dam missing genotype data
    print(paste("WARNING:", sum(tcj$MISS_PARENT_COUNT == 2), "missing elements of the X matrix have been replaced with", miss.x.and.y))
    
    tcj$SAMPLE_ID <- tcj$FAMILY_ID
    # Generate X matrices (assume X matrix is the same for each pool)
    tmp.X.mat <- acast(tcj, SNP_ID ~ SAMPLE_ID , value.var = "B_TRANS_PROB")
    tmp.X.mat <- tmp.X.mat[order(rownames(tmp.X.mat), decreasing = FALSE), ] #order by SNP_ID 
    
    X.mat[[fam.set.combn]] <- tmp.X.mat
    
  }
  
  return(X.mat)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

beta.fun <- function(X.mat, 
                     fkj.and.weight, 
                     fams, 
                     miss.x.and.y = 0.5, 
                     miss.w = 0) {
  #Estimates family contributions to pools based on the method outlined in Henshall et al. 2014
  
  # fkj.and.weight    Data frame. From the output of fkj.and.weight.fun.
  #          1. SNP_ID           is the SNP identifier, 
  #          2. SAMPLE_ID          is the pool identifier
  #          3. FREQ_POOL        is the estimated B allele frequency of the the SNP in the pool (y vector of Henshall et al. 2014 page 6)
  #          4. FREQ_POOL_ERROR_WT is the weight attached to the estimated frequency (w vector of Henshall et al. 2014 page 6)
  
  # miss.x.and.y Number. Value to replace missing values in the X matrix and y vector (Henshall et al. 2014 page 6).  
  
  # miss.w     Number.   Replace missing values in the w vector (Henshall et al. 2014 page 6).  
  
  #Returns##########################################
  
  # beta: Data frame. 
  #          1. SAMPLE_ID   is the pool identifier 
  #          2. SIRE_ID   is the sire identifier, 
  #          3. DAM_ID    is the dam identifier, 
  #          4. FAMILY_ID    is the family identifier. 
  #          5. BETA_STAR is the estimated family contribution (proportion) to the pool as 
  #                       estimated using the pcls function of the mgcv package. Contributions 
  #                       may sum to a number greater than one (see BETA_HAT). See Henshall 
  #                       et al. {, 2014 #2945} page 6.
  #          6. BETA_HAT  is the estimated family contribution (proportion) to the pool 
  #                       where family contributions to each pool are adjusted to sum to one.  
  #                       Calculated as BETA_STAR / (sum of BETA_STAR within each pool).  
  #                       See Henshall et al. {, 2014 #2945} page 6.
  
  print("Running beta.fun")
  
  #Get requried packages
  
 # if("mgcv" %in% installed.packages()[, "Package"] == FALSE) {install.packages("mgcv")} 
  library(mgcv)
  
  fkj.and.weight$SNP_ID           <- as.character(fkj.and.weight$SNP_ID)
  fkj.and.weight$SAMPLE_ID          <- as.integer(fkj.and.weight$SAMPLE_ID)
  fkj.and.weight$FREQ_POOL        <- as.numeric(fkj.and.weight$FREQ_POOL) 
  fkj.and.weight$FREQ_POOL_ERROR_WT <- as.numeric(fkj.and.weight$FREQ_POOL_ERROR_WT) 
  fkj.and.weight <- fkj.and.weight[,c("SNP_ID", "SAMPLE_ID", "FREQ_POOL", "FREQ_POOL_ERROR_WT")]
  
  # Check the list of SNPs the same in input files
  
  if(sum(!unique(fkj.and.weight[,"SNP_ID"]) %in% rownames(X.mat))>0 & 
     sum(!rownames(X.mat) %in% unique(fkj.and.weight[,"SNP_ID"]))>0) {
    print("SNP identifiers do not match in fkj.and.weight and X.mat")
    stop()
  }
  
  # Check that there are more SNP than families
  
  if(nrow(X.mat) <= ncol(X.mat)) {
    print("There are more families than SNP.  Reduce the number of families or increase the number of SNP")
    stop()
  }
  
  #Loop through pools
  fkj.and.weight <- fkj.and.weight[order(fkj.and.weight[,"SNP_ID"], decreasing = FALSE), ]   #order by SNP_ID 
  colnames(fkj.and.weight)   <- c("SNP_ID", "SAMPLE_ID", "FREQ_POOL", "FREQ_POOL_ERROR_WT")  
  
  #Replace missing values
  print(paste("WARNING:", sum(is.na(fkj.and.weight[,"FREQ_POOL"])), "missing elements of the y vector have been replaced with", miss.x.and.y))
  fkj.and.weight[is.na(fkj.and.weight[,"FREQ_POOL"]),"FREQ_POOL"] <- miss.x.and.y
  
  print(paste("WARNING:", sum(is.na(fkj.and.weight[,"FREQ_POOL_ERROR_WT"])), "missing elements of the w vector have been replaced with", miss.w))
  fkj.and.weight[is.na(fkj.and.weight[,"FREQ_POOL_ERROR_WT"]),"FREQ_POOL_ERROR_WT"] <- miss.w
  
  #Empty matrices
  pools <- unique(fkj.and.weight$SAMPLE_ID)
  beta.star   <- matrix(NA, 
                        nrow = length(pools),
                        ncol = ncol(X.mat))
  rownames(beta.star) <- pools
  colnames(beta.star) <- colnames(X.mat)
  
  beta.hat <- beta.star
  
  for (pool in pools) {
    
    beta.star.tmp <- NA
    beta.hat.tmp  <- NA
    
    print(pool)
    
    #Get y and w vectors (see Henshall et al. (2014) pages 5 and 6)
    y <- fkj.and.weight[fkj.and.weight[,"SAMPLE_ID"] == pool, "FREQ_POOL"]
    w <- matrix(fkj.and.weight[fkj.and.weight[,"SAMPLE_ID"] == pool, "FREQ_POOL_ERROR_WT"],ncol = 1)
    
    # Estimate beta using pcls function of mgcv package
    
    # Inputs for pcls funtion. See help(pcls)
    M <- list( y   = y,
               w   = w,
               X   = X.mat,
               C   = matrix(0,0,0),
               p   = rep(1, ncol(X.mat)),
               off = array(0,0),
               S   = list(),
               sp  = array(0,0),
               Ain = diag(ncol(X.mat)),
               bin = rep(0, ncol(X.mat)) )
    
    
    try( beta.star.tmp     <- pcls(M) )
    beta.star.tmp <- round(beta.star.tmp,10)
    
    beta.hat.tmp      <- beta.star.tmp / sum(beta.star.tmp)
    
    beta.star[rownames(beta.star) == pool,] <- beta.star.tmp
    beta.hat[rownames(beta.hat) == pool,]   <- beta.hat.tmp
    
  }
  
  #Generate beta data frame
  
  #Rearrange beta.star
  tmp.beta.star <- NULL
  #Loop through columns
  for (c in 1:ncol(beta.star)) {
    tmp.1 <- data.frame(SAMPLE_ID = rownames(beta.star),
                        FAMILY_ID = colnames(beta.star)[c],
                        BETA_STAR = beta.star[,c])
    tmp.beta.star <- rbind(tmp.beta.star,tmp.1)
  }
  tmp.beta.star[,"SAMPLE_ID"] <- as.integer(as.character(tmp.beta.star[,"SAMPLE_ID"]))
  tmp.beta.star[,"FAMILY_ID"] <- as.integer(as.character(tmp.beta.star[,"FAMILY_ID"]))
  
  #Rearrange beta.hat
  tmp.beta.hat <- NULL
  #Loop through columns
  for (c in 1:ncol(beta.hat)) {
    tmp.1 <- data.frame(SAMPLE_ID = rownames(beta.hat),
                        FAMILY_ID = colnames(beta.hat)[c],
                        BETA_HAT = beta.hat[,c])
    tmp.beta.hat <- rbind(tmp.beta.hat,tmp.1)
  }
  tmp.beta.hat[,"SAMPLE_ID"] <- as.integer(as.character(tmp.beta.hat[,"SAMPLE_ID"]))
  tmp.beta.hat[,"FAMILY_ID"]  <- as.integer(as.character(tmp.beta.hat[,"FAMILY_ID"]))
  
  beta <- merge(tmp.beta.star, tmp.beta.hat, by = c("SAMPLE_ID", "FAMILY_ID"))
  
  #  beta <- merge(beta, fams, by = "FAMILY_ID", all.x = TRUE)
  beta <- left_join(beta, fams, by = "FAMILY_ID")
  
  beta <- beta[,c("SAMPLE_ID", "SIRE_ID", "DAM_ID", "FAMILY_ID", "BETA_STAR", "BETA_HAT")] 
  
  beta <- beta[order(beta$SAMPLE_ID, decreasing = FALSE), ] 
  
  return(beta)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

beta.min.ss.fun <- function(X.mat, fkj.and.weight, miss.x.and.y = 0.5, miss.w = 0, fam.set.combns, fams) {
  #Estimates family contributions to pools based on the method outlined in Henshall et al. 2014.  And then
  #estimates beta constrained to possible values given the known number of individuals in a pool an minimum sum of squares identified
  #and possible contributors.  That is elements of BETA_MIN_SS must be 0 or a multiple of 1 / (contributors to the pool)
  
  # fkj.and.weight    Data frame. From the output of fkj.and.weight.fun.
  #          1. SNP_ID           is the SNP identifier, 
  #          2. SAMPLE_ID          is the pool identifier
  #          3. FREQ_POOL        is the estimated B allele frequency of the the SNP in the pool (y vector of Henshall et al. 2014 page 6)
  #          4. FREQ_POOL_ERROR_WT is the weight attached to the estimated frequency (w vector of Henshall et al. 2014 page 6)
  
  # miss.x.and.y Number. Value to replace missing values in the X matrix and y vector (Henshall et al. 2014 page 6).  
  
  # miss.w     Number.   Replace missing values in the w vector (Henshall et al. 2014 page 6).  
  
  # n.by.pool  Data frame. See indiv.ped.file.fun. Generally NULL but if not:
  #          1. SAMPLE_ID  is the pool identifier and 
  #          2. N_INDIV is the number of individuals that contributed to the pool.    
  
  #Returns##########################################
  
  # beta: Data frame. 
  #          1. SAMPLE_ID   is the pool identifier 
  #          2. SIRE_ID   is the sire identifier, 
  #          3. DAM_ID    is the dam identifier, 
  #          4. FAMILY_ID    is the family identifier. 
  #          5. BETA_STAR is the estimated family contribution (proportion) to the pool as 
  #                       estimated using the pcls function of the mgcv package. Contributions 
  #                       may sum to a number greater than one (see BETA_HAT). See Henshall 
  #                       et al. {, 2014 #2945} page 6.
  #          6. BETA_HAT  is the estimated family contribution (proportion) to the pool 
  #                       where family contributions to each pool are adjusted to sum to one.  
  #                       Calculated as BETA_STAR / (sum of BETA_STAR within each pool).  
  #                       See Henshall et al. {, 2014 #2945} page 6.
  #          7. BETA_MIN_SS beta constrained to possible contbributions given n.by.pool and min sum of squares identified.  That is
  #                       elements of BETA_MIN_SS must be 0 or a multiple of 1 / (contributors to the pool)
  
  print("Running beta.min.ss.fun")
  
  #Get requried packages
  
 # if("mgcv" %in% installed.packages()[, "Package"] == FALSE) {install.packages("mgcv")} 
  library(mgcv)
  
  fkj.and.weight$SNP_ID             <- as.character(fkj.and.weight$SNP_ID)
  fkj.and.weight$SAMPLE_ID          <- as.integer(fkj.and.weight$SAMPLE_ID)
  fkj.and.weight$FREQ_POOL          <- as.numeric(fkj.and.weight$FREQ_POOL) 
  fkj.and.weight$FREQ_POOL_ERROR_WT <- as.numeric(fkj.and.weight$FREQ_POOL_ERROR_WT) 
  fkj.and.weight <- fkj.and.weight[,c("SNP_ID", "SAMPLE_ID", "FREQ_POOL", "FREQ_POOL_ERROR_WT")]
  
  # Check the list of SNPs the same in input files
  
  if(sum(!unique(fkj.and.weight[,"SNP_ID"]) %in% rownames(X.mat))>0 & 
     sum(!rownames(X.mat) %in% unique(fkj.and.weight[,"SNP_ID"]))>0) {
    print("SNP identifiers do not match in fkj.and.weight and X.mat")
    stop()
  }
  
  # Check that there are more SNP than families
  if(nrow(X.mat) <= ncol(X.mat)) {
    print("There are more families than SNP.  You need to reduce the number of families or increase the number of SNP to run Least_squares method")
    stop()
  }
  
  #Loop through pools
  fkj.and.weight <- fkj.and.weight[order(fkj.and.weight[,"SNP_ID"], decreasing = FALSE), ]   #order by SNP_ID 
  colnames(fkj.and.weight)   <- c("SNP_ID", "SAMPLE_ID", "FREQ_POOL", "FREQ_POOL_ERROR_WT")  
  
  #Replace missing values
  print(paste("WARNING:", sum(is.na(fkj.and.weight[,"FREQ_POOL"])), "missing elements of the y vector have been replaced with", miss.x.and.y))
  fkj.and.weight[is.na(fkj.and.weight[,"FREQ_POOL"]),"FREQ_POOL"] <- miss.x.and.y
  
  print(paste("WARNING:", sum(is.na(fkj.and.weight[,"FREQ_POOL_ERROR_WT"])), "missing elements of the w vector have been replaced with", miss.w))
  fkj.and.weight[is.na(fkj.and.weight[,"FREQ_POOL_ERROR_WT"]),"FREQ_POOL_ERROR_WT"] <- miss.w
  
  #Empty matrices
  pools <- unique(fkj.and.weight$SAMPLE_ID)
  beta.star   <- matrix(NA, 
                        nrow = length(pools),
                        ncol = ncol(X.mat))
  rownames(beta.star) <- pools
  colnames(beta.star) <- colnames(X.mat)
  
  beta.hat <- beta.star
  beta.contstrained <- beta.star
  beta.contstrained <- as.data.frame(beta.contstrained)
  
  for (pool in pools) {
    
    beta.star.tmp <- NA
    beta.hat.tmp  <- NA
    
    print(paste("Pool",pool))
    
    #Get y and w vectors (see Henshall et al. (2014) pages 5 and 6)
    y <- fkj.and.weight[fkj.and.weight[,"SAMPLE_ID"] == pool, "FREQ_POOL"]
    w <- matrix(fkj.and.weight[fkj.and.weight[,"SAMPLE_ID"] == pool, "FREQ_POOL_ERROR_WT"],ncol = 1)
    
    # Estimate beta using pcls function of mgcv package
    
    # Inputs for pcls funtion. See help(pcls)
    M <- list( y   = y,
               w   = w,
               X   = X.mat,
               C   = matrix(0,0,0),
               p   = rep(1, ncol(X.mat)),
               off = array(0,0),
               S   = list(),
               sp  = array(0,0),
               Ain = diag(ncol(X.mat)),
               bin = rep(0, ncol(X.mat)) )
    
    try( beta.star.tmp     <- pcls(M) )
    beta.star.tmp <- round(beta.star.tmp,10)
    
    beta.hat.tmp      <- beta.star.tmp / sum(beta.star.tmp)
    
    beta.star[rownames(beta.star) == pool,] <- beta.star.tmp
    beta.hat[rownames(beta.hat) == pool,]   <- beta.hat.tmp
    
    #compute beta.contstrained
    
    #Get number of contributors in the pool
    # n <- n.by.pool[n.by.pool[,"SAMPLE_ID"] == pool,"N_INDIV"]
    n.in.pool <- length(unique(fam.set.combns[,"FAM_SET_ID"]))
    
    #estimate integer contributions of families in pool
    #    beta.integer <- round(beta.hat.tmp*n)
    #    
    #get all possible combinations of beta constrained to beta.integer +/- tmp
    #    tmp <- max(2, round(0.15*n))
    #    beta.integer <- data.frame(MIN = beta.integer - 2,
    #                               EXPECTED = beta.integer,
    #                               MAX = beta.integer + 2)
    #    beta.integer[beta.integer[,] < 0] <- 0
    #    beta.integer[beta.integer[,] > n] <- n
    
    #    combinations <- NULL
    #    for (i in 1:nrow(beta.integer)) {
    #      if(i == 1) {
    #        combinations <- as.data.frame(beta.integer[i,"MIN"]:beta.integer[i,"MAX"])
    #      } else {
    #        combinations <- merge(combinations,(beta.integer[i,"MIN"]:beta.integer[i,"MAX"]))
    #        combinations <- combinations[rowSums(combinations) <= n,] #remove row if row sum is greater than n
    #        colnames(combinations) <- 1:ncol(combinations)
    #      }
    #    }
    
    #remove rows that don't sum to n
    #    combinations <- as.matrix(combinations[rowSums(combinations) == n,])
    #    combinations <- combinations / n #express as proportion
    
    combinations <- parent.combns.fun(fams, 
                                      n.in.pools = n.in.pool,
                                      fam.set.combns = fam.set.combns)$parent.combns
    combinations <- combinations[,1:(n.in.pool+1)]
    
    #compute ||(sqrt(w) * (X * beta - y))||^2
    
    combinations[,"SUM_OF_SQUARES"] <- NA
    for (j in 1:nrow(combinations)) {
      print(paste("Make beta.min.ss = FALSE if this is too slow.  Sum of squares", j, "of", nrow(combinations)))
      tmp.beta <- combinations[j,!colnames(combinations) %in% c("FAM_COMBN_ID", "SUM_OF_SQUARES")]
      tmp.beta <- colnames(X.mat) %in% tmp.beta
      tmp.beta <- tmp.beta / n.in.pool
      # tmp.beta <- as.matrix(combinations[j,-ncol(combinations)])
      combinations[j,"SUM_OF_SQUARES"] <- sum((sqrt(w) * (X.mat %*% tmp.beta - y))^2)
    }
    
    #identify minimum SUM_OF_SQUARES
    beta.contstrained.tmp <- combinations[combinations[,"SUM_OF_SQUARES"] == min(combinations[,"SUM_OF_SQUARES"]),-ncol(combinations)]
    beta.contstrained.tmp <- colSums(beta.contstrained.tmp) / nrow(beta.contstrained.tmp) #in case min SUM_OF_SQUARES is the same for multiple combinations
    
    beta.contstrained.tmp <- colnames(beta.contstrained) %in% beta.contstrained.tmp[-1] / n.in.pool
    
    beta.contstrained[rownames(beta.contstrained) == pool,]   <- as.vector(beta.contstrained.tmp)
  }
  
  #Generate beta data frame
  
  #Rearrange beta.star
  tmp.beta.star <- NULL
  #Loop through columns
  for (c in 1:ncol(beta.star)) {
    tmp.1 <- data.frame(SAMPLE_ID = rownames(beta.star),
                        FAMILY_ID = colnames(beta.star)[c],
                        BETA_STAR = beta.star[,c])
    tmp.beta.star <- rbind(tmp.beta.star,tmp.1)
  }
  tmp.beta.star[,"SAMPLE_ID"] <- as.integer(as.character(tmp.beta.star[,"SAMPLE_ID"]))
  tmp.beta.star[,"FAMILY_ID"] <- as.integer(as.character(tmp.beta.star[,"FAMILY_ID"]))
  
  #Rearrange beta.hat
  tmp.beta.hat <- NULL
  #Loop through columns
  for (c in 1:ncol(beta.hat)) {
    tmp.1 <- data.frame(SAMPLE_ID = rownames(beta.hat),
                        FAMILY_ID = colnames(beta.hat)[c],
                        BETA_HAT = beta.hat[,c])
    tmp.beta.hat <- rbind(tmp.beta.hat,tmp.1)
  }
  tmp.beta.hat[,"SAMPLE_ID"] <- as.integer(as.character(tmp.beta.hat[,"SAMPLE_ID"]))
  tmp.beta.hat[,"FAMILY_ID"]  <- as.integer(as.character(tmp.beta.hat[,"FAMILY_ID"]))
  
  
  #Rearrange beta.contstrained
  tmp.beta.contstrained <- NULL
  #Loop through columns
  for (c in 1:ncol(beta.contstrained)) {
    tmp.1 <- data.frame(SAMPLE_ID = rownames(beta.contstrained),
                        FAMILY_ID = colnames(beta.contstrained)[c],
                        BETA_MIN_SS = beta.contstrained[,c])
    tmp.beta.contstrained <- rbind(tmp.beta.contstrained,tmp.1)
  }
  tmp.beta.contstrained[,"SAMPLE_ID"] <- as.integer(as.character(tmp.beta.contstrained[,"SAMPLE_ID"]))
  tmp.beta.contstrained[,"FAMILY_ID"] <- as.integer(as.character(tmp.beta.contstrained[,"FAMILY_ID"]))
  
  
  beta <- merge(tmp.beta.star, tmp.beta.hat, by = c("SAMPLE_ID", "FAMILY_ID"))
  beta <- merge(beta, tmp.beta.contstrained, by = c("SAMPLE_ID", "FAMILY_ID"))
  
  #beta <- merge(beta, fams, by = "FAMILY_ID", all.x = TRUE)
  beta <- left_join(beta, fams, by = "FAMILY_ID")
  
  beta <- beta[,c("SAMPLE_ID", "SIRE_ID", "DAM_ID", "FAMILY_ID", "BETA_STAR", "BETA_HAT", "BETA_MIN_SS")] 
  
  beta <- beta[order(beta$SAMPLE_ID, decreasing = FALSE), ] 
  
  return(beta)
}

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

#' @export

bar.plot.fun <- function(beta,
                         file.name = "",
                         var = "BETA_HAT",
                         heading = "Estimated family contributions to pooled samples",
                         plot.to.heading.height = 20,
                         font.size.heading = 3,
                         font.size.y.axis = 2,
                         font.size.x.axis = 2 
) {
  
  #Generates abar plot of contributions of families to pools and dendrogram that groups pools
  
  #Args##########################################
  # beta: Data frame.  From output of beta.fun. 
  #          1. SAMPLE_ID   is the pool identifier 
  #          2. FAMILY_ID    is the family identifier. 
  #          3. BETA_STAR is the estimated family contribution (proportion) to the pool as 
  #                       estimated using the pcls function of the mgcv package. Contributions 
  #                       may sum to a number greater than one (see BETA_HAT). See Henshall 
  #                       et al. {, 2014 #2945} page 6.
  #          4. BETA_HAT  is the estimated family contribution (proportion) to the pool 
  #                       where family contributions to each pool are adjusted to sum to one.  
  #                       Calculated as BETA_STAR / (sum of BETA_STAR within each pool).  
  #                       See Henshall et al. {, 2014 #2945} page 6.
  
  # file.name:              Text. Name of thebar plot file
  # var:                    Text.  Variable to plot "BETA_STAR" or "BETA_HAT"
  # heading:                Text. Title of thebar plot 
  # plot.to.heading.height: Number. Height of the title relative to the height of thebar plot 
  # font.size.heading:      Number. Font size ofbar plot heading.
  # font.size.y.axis:       Number. Font size ofbar plot y axis labels
  # font.size.x.axis:       Number. Font size ofbar plot x axis labels
  
  # Returns:
  # heat.map.png saved to the working directory
  
  print("Running bar.plot.fun")
  
  #Load required packages
 # if("RColorBrewer" %in% installed.packages()[, "Package"] == FALSE) {install.packages("RColorBrewer")} 
  library(RColorBrewer)
  
 # if("gplots" %in% installed.packages()[, "Package"] == FALSE) {install.packages("gplots")} 
  library(gplots)
  
 # if("reshape2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("reshape2")} 
  library(reshape2)  
  
 # if("ggplot2" %in% installed.packages()[, "Package"] == FALSE) {install.packages("ggplot2")} 
  library(ggplot2)
  
  #Name columns and assign class
  beta$SAMPLE_ID     <- as.integer(beta$SAMPLE_ID  )
  beta$FAMILY_ID      <- as.integer(beta$FAMILY_ID   )  
  beta$BETA_STAR    <- as.numeric(beta$BETA_STAR )
  beta$BETA_HAT    <- as.numeric(beta$BETA_HAT )
  
  if("BETA_MIN_SS" %in% colnames(beta)) {
    if(sum(is.na(beta[,"BETA_MIN_SS"])) != nrow(beta)) { #if BETA_MIN_SS all equal NA
      beta$BETA_MIN_SS    <- as.numeric(beta$BETA_MIN_SS)
      beta$MIN_SUM_SQ <- beta$BETA_MIN_SS != 0
    } else {
      beta <- beta[,c("SAMPLE_ID", "FAMILY_ID", "BETA_STAR", "BETA_HAT")]
      beta$MIN_SUM_SQ <- FALSE
    }
  } else {
    beta <- beta[,c("SAMPLE_ID", "FAMILY_ID", "BETA_STAR", "BETA_HAT")]
    beta$MIN_SUM_SQ <- FALSE
  }
  
  beta$VAR <- beta[,colnames(beta) == var]
  beta$FAMILY_ID <- as.character(beta$FAMILY_ID)
  
  # Bar graph, family on x-axis, color fill grouped by pool -- use position_dodge()
  barplot <- ggplot(data=beta, aes(x=FAMILY_ID, y=VAR, fill=MIN_SUM_SQ)) 
  #barplot <- barplot + ylab(var) 
  barplot <- barplot + labs(y = var)
  barplot <- barplot + geom_bar(stat="identity", position=position_dodge())
  barplot <- barplot + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) 
  ggsave(filename = paste(file.name,".bar.png",sep=""), plot = barplot)
  
}
