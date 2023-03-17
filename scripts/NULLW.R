###################################################################################################################################
################################################repeated association test with null-w##############################################
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
# argument5: top candidate list from TOPCAN script
# argument6: 0.05,  the threshold for q-value in null-W results (Empirical p-values of null-W were converted to q-values for estimation of False Discovery Rate)
# argument7: save directory
# argument8: Number of threads for parallel computation 

# example: Rscript NULLW.R annuus_GWAS_result.txt argophyllus_GWAS_result.txt P 5kbwindow_recombination 10 0.9999 ~ 150
TEST <- toupper(args[3])
recombination_table <- args[4]
TOPCANS <- args[5]
q <- as.numeric(args[6])
save_dir <- args[7]
nthreads = as.numeric(args[8])

################################################################Output#################################################################################################

#$species1_$species2_$variablename_nullw_result #null-W test results for species 2 topcandidates in species 1 genome
#$species1_$species2_$variablename_nullw_result_FDR_$q #null-W test results for species 2 topcandidates in species 1 genome which passed the FDR test
#$species2_$species_1$variablename_nullw_result #null-W test results for species 1 topcandidates in species 2 genome
#$species2_$species1_$variablename_nullw_result_FDR_$q #null-W test results for species 1 topcandidates in species 2 genome which passed the FDR test

#######################################################################data preperation for null-w test################################################################
registerDoParallel(cores=nthreads)

fread(TOPCANS,
      showProgress = F,
      header = T) -> data #top candidate analysis output

data %>%
  distinct(species) %>%
  pull() -> species_pair #species pair to work on

rbind(data %>%
        group_by(window) %>%
        tally() %>%
        filter(n == 2) %>%
        select(window) %>%
        mutate(species = species_pair[1]) %>%
        mutate(species2 = species_pair[2]),
      data %>%
        group_by(window) %>%
        tally() %>%
        filter(n == 2) %>%
        select(window) %>%
        mutate(species = species_pair[2]) %>%
        mutate(species2 = species_pair[1])) -> orth_window_species # list of otholog windows


data %>% 
  select(species2 = species ,
         variable ,
         window ,
         species2_top_candidate = top_candidate) -> top_compare#list of all windows in both species #331219

ortho_species1.2_top <-foreach(i=1:2, .combine='rbind', .errorhandling='stop') %dopar% {
  fread(args[i],
        showProgress = F,
        header = F) %>% #3332702
    select(species = V1 ,
           variable = V2 ,
           snp_ID = V3 ,
           window = V4 ,
           Value = V5) %>%
    left_join(.,
              data) %>%
    filter(!is.na(snp_count)) %>% #add top candidate and SNP_count window info to each and window #3332702
    select(-outlier_count ,
           -Qbinom) %>%
    left_join(.,
              orth_window_species) %>% #add ortholog daa
    filter(!is.na(species2)) %>% #filter for non ortholog window #2990417
    left_join(.,
              top_compare) %>% #data to show if the ortholog window is top candidate in the other species
    filter(!is.na(species2_top_candidate)) 
}

rm(data,orth_window_species,top_compare)


for (i in c(1,2)) { #species
  
  ortho_species1.2_top %>%
    filter(species == species_pair[i]) %>%
    select(window) %>%
    distinct(window) %>%
    left_join(.,
              rename(
                fread(recombination_table ,
                      sep="\t" ,
                      header = F) , 
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
                                 include.lowest=TRUE)) -> window_recom_ortho # define recombination bins
  ortho_species1.2_top %>%
    filter(species == species_pair[i]) %>%
    select(variable,Value,snp_ID,window,snp_count,top_candidate,species2_top_candidate) %>% 
    left_join(.,window_recom_ortho) %>%
    filter(!is.na(recomb_rate)) %>%
    select(-recomb_rate) -> two_variable_orth #assign recombination bins to each window
  
  rm(window_recom_ortho)
  
  two_variable_orth %>%
    distinct(recomb_quantile) %>%
    pull(recomb_quantile) -> quantile_list # list of recombination bins
  

  Final_result <- NULL
  
  for (z in 1:length(quantile_list)) { #recombination bins
    
    top_list <- two_variable_orth %>% 
      filter(recomb_quantile == quantile_list[z]) %>%
      filter(top_candidate=="top") %>% 
      distinct(window) %>%
      pull(window) #top candidate windows for species 1
    
    background_10_ID <- two_variable_orth %>%
      filter(recomb_quantile == quantile_list[z]) %>%
      filter(!window %in% top_list) %>%
      distinct(snp_ID) %>%
      sample_n(10000) %>%
      pull(snp_ID) #1000 random snp IDs in non-top candidate windows for species 1
    
    background_10KB <- two_variable_orth %>%
      filter(recomb_quantile == quantile_list[z]) %>%
      filter(snp_ID %in% background_10_ID) %>%
      pull(Value) #values for 1000 random snp in non-top candidate windows for species 1
    
    window_list <- two_variable_orth %>%
      filter(recomb_quantile == quantile_list[z]) %>%
      filter(!window %in% top_list) %>%
      filter(!is.na(Value)) %>%
      filter(snp_count > 5) %>%
      distinct(window) %>%
      pull(window) #non top candidate windows with more than 5 SNPs in species 1
    
    Values <- two_variable_orth %>%
      filter(recomb_quantile == quantile_list[z]) %>%
      select(window,Value) # all values for windows in species 1
    
    window_list1 <- two_variable_orth %>%
      filter(recomb_quantile == quantile_list[z]) %>%
      filter(species2_top_candidate=="top") %>%
      filter(snp_count > 3) %>%
      distinct(window) %>%
      pull(window) #windows for species 1 with more than 3 snps that are top candidate in species2 
    
    N_window <- two_variable_orth %>%
      filter(recomb_quantile == quantile_list[z]) %>%
      filter(species2_top_candidate == "top") %>%
      distinct(window) %>%
      nrow() #number of top candidate windows for species 2
    
    
    if (TEST == "P") {
      
      background_10KB %>%
        log10(.)*-1 -> background_10KB
      
      Values %>%
        mutate(Value=log10(Value)*-1) -> Values
      
      print("P values are adjusted")
      
    } else if (TEST == "R") { 
      
      background_10KB^2 -> background_10KB
      
      Values %>%
        mutate(Value=Value^2) -> Values
      print("R values are adjusted")
      
    } else if (TEST == "BF") { 
      print("BF values are adjusted")
    } else { 
      print("Wrong value in the third argument")
    }
    
    
    if (length(window_list1) > 0) { 
      
      nulldist<-foreach(k = 1:length(window_list), .combine='rbind', .errorhandling='stop') %dopar% {#null distribution
        
        data.frame(N= sum(is.na(pull(filter(Values,window==window_list[k]),Value))==F),
                   W=as.vector(wilcox.test(background_10KB,pull(filter(Values,window==window_list[k]),Value))$statistic)
                   ) %>% 
          filter(!is.na(W)) %>%                                                                                           
          mutate(Zscore= (2*W-10000*N) / sqrt(10000*N*(10000+N+1)/3)) %>%
          select(Zscore)
      }#null distribution
      
      rm(window_list)

      Result <- foreach(l=1:length(window_list1), .combine='rbind', .errorhandling='stop') %dopar% {
        
        data.frame(species1=species_pair[i],
                   species2=species_pair[species_pair != species_pair[i]],
                   Variable=distinct(two_variable_orth,variable),
                   recomb_quantile=as.character(quantile_list[z]),
                   Window=window_list1[l],
                   W=as.vector(wilcox.test (pull(filter(Values,window==window_list1[l]),Value),background_10KB)$statistic),#wilcox test, x=(values of first not top window for the second species), y=(10000 p_values in first variable that are not in top windows for first species)
                   P=wilcox.test (pull(filter(Values,window==window_list1[l]),Value),background_10KB)$p.value,
                   Mean_test=mean(pull(filter(Values,window==window_list1[l]),Value)),
                   N_samples=length(pull(filter(Values,window==window_list1[l]),Value)[!is.na(pull(filter(Values,window==window_list1[l]),Value))]),
                   Mean_BG=mean(background_10KB),
                   N_BG=length(background_10KB[!is.na(background_10KB)]),
                   Count=N_window)
      }
      
      temporary_result <- Result %>% mutate(Pred_P=(2* W - N_BG * N_samples)/sqrt(N_BG * N_samples*(N_BG+N_samples+1)/3))
      
      rm(Values,window_list1,background_10KB)
      
      Emp_result <- foreach(m=1:nrow(temporary_result), .combine='rbind', .errorhandling='stop') %dopar% {
        data.frame(Emp_p = 1-(sum(as.numeric(temporary_result[m,13]) > pull(nulldist,Zscore))/nrow(nulldist)),Window=temporary_result[m,5])
      }
      
      left_join(temporary_result
                ,Emp_result) %>%
        filter(!is.na(Emp_p)) %>%
        rbind(. ,
              Final_result) -> Final_result
      
      rm(nulldist,Result,temporary_result,Emp_result,N_window)} else {
        data.frame(variable = distinct(two_variable_orth,variable) ,
                   recomb_quantile = quantile_list[z] ,
                   species1 = species_pair[i] ,
                   species2 = species_pair[species_pair != species_pair[i]]) %>%
          fwrite(file = paste(save_dir,
                              "/",
                              species_pair[1] ,
                              "_",
                              species_pair[2] ,
                              "_",
                              distinct(two_variable_orth,variable) ,
                              "_nullw_result_No_outlier") ,
                 append = T ,
                 sep = "\t")
        rm(background_10KB,window_list,Values,window_list1,N_window)}
    
  }#quantile loop 

#################################################################################q-value calculation and FDR test################################################################

 Final_result$q_value <- p.adjust(Final_result$Emp_p, method = "BH")

fwrite(Final_result ,
       file = paste0(save_dir,
                    "/",
                    species_pair[i] ,
                    "_",
                    species_pair[species_pair != species_pair[i]] ,
                    "_",
                    distinct(two_variable_orth,variable) ,
                    "_nullw_result") ,
       append = T ,
       sep = "\t")

fwrite(filter(Final_result,q_value <= q) ,
       file = paste0(save_dir,
                    "/",
                    species_pair[i] ,
                    "_",
                    species_pair[species_pair != species_pair[i]] ,
                    "_",
                    distinct(two_variable_orth,variable) ,
                    "_nullw_result_FDR" ,
                    "_",
                    q ) ,
       append = T ,
       sep = "\t")

rm(Final_result)
}



