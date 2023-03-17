###################################################################################################################################
################################################Identification of top candidate windows ###########################################
###################################################################################################################################


###############################################################Packages############################################################
library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
#############################################################Arguments##############################################################
args = commandArgs(trailingOnly = TRUE)

# argument1: association result (GWAS,BYPASS, or correlation) in species 1, includes for columns (with out header) as follow: 1. species  2. Variable 3.SNP_ID (example Ha412HOChr01:180057) 4.window ID (example Ha412HOChr01:180001-185000) 5.asociation value (P,BF or r)
# argument2: association result (GWAS,BYPASS, or correlation) in species 2, includes for columns (with out header) as follow: 1. species  2. Variable 3.SNP_ID (example Ha412HOChr01:180057) 4.window ID (example Ha412HOChr01:180001-185000) 5.asociation value (P,BF or r)
# argument3: association value , P for P value, R for correlation, BF for Bayes Factor
# argument4: recombination ratio file "5kbwindow_recombination"

# argument5:  #The threshold for the significant association. 
#Example: for P-value the 0.01 cut off means the bottom one percent quantile of the association tests (smaller P-value = stronger association)
#for correlation, the threshold should be 0.99 (as larger correlation means stronger association), 
#for BF just enter the threshold for example 10, it will filter for values larger than 10

# argument6: #The threshold for identifying top candidate windows (0.9999 quantile of the binomial expectation)
# argument7: save directory
# argument8: Number of threads for parallel computation 

# example: Rscript TOPCAN.R annuus_GWAS_result.txt argophyllus_GWAS_result.txt P 5kbwindow_recombination 10 0.9999 ~ 150
TEST <- toupper(args[3])
recombination_table <- args[4]
association_cut <- as.numeric(args[5])
dbinom_cut <- as.numeric(args[6])
save_dir <- args[7]
nthreads = as.numeric(args[8])

################################################################Output##############################################################

#$save_dir_$argument1_topcandidate #the file indicates if each window is a top candidate for the species in argument 1, and contains information on number of SNPs (snp_count) and number of significantly associated snps(outlier_count) for each window
#$save_dir_$argument2_topcandidate #the file indicates if each window is a top candidate for the species in argument 2, and contains information on number of SNPs (snp_count) and number of significantly associated snps(outlier_count) for each window
#$save_dir_$argument1_$argument2_topcandidate #concatenated version for both spices

#################################################Identification of top candidate windows######################################################################################

registerDoParallel(cores=nthreads)
data <- foreach(i=1:2, .combine='rbind', .errorhandling='stop') %dopar% {
  fread(args[i],
        showProgress = F,
        header = F) %>%
    mutate(test= TEST) %>%
    mutate(value = ifelse(test == "R",
                           V5^2,
                           V5)) %>%
    select(species = V1,
           variable = V2,
           window = V4,
           value) %>%
    mutate(outlier = ifelse(value >= ifelse(TEST == "BF", association_cut, quantile(as.numeric(value),association_cut)),
                            "outlier",
                            "non_outlier")) %>% 
    group_by(species,
             variable,
             window,
             outlier) %>% 
    summarise(snp_count = n()) %>%
    ungroup() %>%
    spread(outlier,
           snp_count) %>%
    mutate(outlier = replace_na(outlier, 0) ,
           non_outlier = replace_na(non_outlier,
                                    0)) %>%
    mutate(snp_count = as.numeric(outlier) + as.numeric(non_outlier)) %>%
    select(-non_outlier) %>%
    rename(outlier_count = outlier) %>%
    left_join(. , 
              rename(fread(recombination_table,sep="\t", header = F) , 
                     window = V1,
                     recomb_rate = V2)) %>%
    filter(!is.na(recomb_rate)) %>%
    mutate(recomb_quantile = cut(recomb_rate, 
                                 breaks = quantile(recomb_rate ,
                                                   probs=seq(0 ,
                                                             1 , 
                                                             by=0.2) ,
                                                   na.rm=TRUE) , 
                                 labels = c("0-20" ,
                                            "20-40" ,
                                            "40-60" ,
                                            "60-80" ,
                                            "80-100"), 
                                 include.lowest=TRUE)) %>%
    select(-recomb_rate) -> SNP_count
  
  left_join(SNP_count ,
            SNP_count %>% 
              group_by(species ,
                       variable ,
                       recomb_quantile) %>% 
              summarise(totsnp = sum(snp_count) ,
                        totout = sum(outlier_count)) %>% 
              mutate(expected = totout/totsnp) %>% 
              select(-totsnp,
                     -totout)) %>% 
    mutate(Qbinom = qbinom (dbinom_cut ,
                            snp_count ,
                            expected)) %>%
    mutate(top_candidate = ifelse(outlier_count > Qbinom,
                                  "top",
                                  "not")) %>%
    select(species,
           variable,
           window,
           outlier_count,
           snp_count,
           Qbinom,
           top_candidate) %>%
    fwrite(paste0(save_dir,
                  "/",
                  basename(args[i]),
                  "_topcandidate"),
           sep = "\t",
           col.names = T)
  
  rm(SNP_count)
  
  fread(paste0(save_dir,
               "/",
               basename(args[i]),
               "_topcandidate"),
        showProgress = F,
        header = T) 
}


data %>% 
  fwrite(paste0(save_dir,
                "/",
                basename(args[1]),
                "_",
                basename(args[2]),
                "_topcandidate"),
         sep = "\t",
         col.names = T)
