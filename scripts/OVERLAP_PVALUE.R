library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(IRanges)

##################################################################Arguments#########################################################################################
args = commandArgs(trailingOnly = TRUE)
HAPLOBLOCKS = args[1] 
FILE = args[2] 
COMPARISON = args[3] 
GENOME = args[4]
SAVE_DIR= args[5]
nthreads = as.numeric(args[6])
#####################################################################################################################################################################
PERMUT = 10000 #number of permutation for calculation of P valuels 


fread(HAPLOBLOCKS,
      header = T) -> haploblocks.bed


Taxa <- unlist(strsplit(COMPARISON, "_"))

registerDoParallel(cores=nthreads)
for (i in 1:length(Taxa)) {
  haploblocks.bed %>%
    filter(species == Taxa[i]) %>%
    select(chr,start,end,haploblock) %>%
    fwrite(paste0(SAVE_DIR,
                  "/",
                  Taxa[i],
                  "_haploblocks.bed"),
           sep = "\t",
           col.names = F,
           quote = F)
  
  foreach(x=1:PERMUT, .combine='rbind', .errorhandling='stop') %dopar% {
    system(paste0("bedtools shuffle -i ",
                  SAVE_DIR,
                  "/",
                  Taxa[i],
                  "_haploblocks.bed",
                  " -g ",
                  GENOME,
                  " -noOverlapping > ",
                  SAVE_DIR,
                  "/",
                  Taxa[i],
                  "_perm",
                  x,
                  "_haploblocks.bed")) 
    
    fread(
      paste0(
        SAVE_DIR,
        "/", 
        Taxa[i],
        "_perm",
        x,
        "_haploblocks.bed"),
      sep = "\t",
      header = F) %>%
      mutate(species = Taxa[i],
             permutation = x) %>%
      select(chr=V1,start=V2,end=V3,species,permutation,sv_name_new=V4) %>%
      fwrite(
        paste0(
          SAVE_DIR,
          "/randomized_haploblocks.bed"),
        sep = "\t",
        col.names = F,
        append = T)
    
    system(paste0("rm ", 
                  SAVE_DIR,
                  "/",
                  Taxa[i],
                  "_perm",
                  x,
                  "_haploblocks.bed"))
  }
  
  system(paste0("rm ", 
                SAVE_DIR,
                "/",
                Taxa[i],
                "_haploblocks.bed"))
}


  fread(
    paste0(
      SAVE_DIR,
      "/randomized_haploblocks.bed")) %>%
    select(chromosome = V1,
           start = V2,
           end = V3,
           species = V4,
           permutation = V5,
           haploblock = V6) -> randomized_haploblocks
  
  system(paste0("rm ", 
                SAVE_DIR,
                "/randomized_haploblocks.bed"))


  
  
  fread(FILE,sep = "\t",header = T) %>%
    select(comparison,
           CRA) %>%
    distinct(CRA,.keep_all = T) %>%
    separate(comparison, into = c("Taxa1","Taxa2"),sep = "_",remove = F) %>%
    separate(CRA,into = c("chromosome","start","end"),sep = ":",remove = F) %>%
    mutate(start=as.numeric(start),end=as.numeric(end)) %>%
    select(chromosome,
           Taxa1,
           Taxa2,
           start,
           end,
           comparison)-> convergent_cluster
  
  perm_result <- NULL
  

    for (x in 1:PERMUT) {
      randomized_haploblocks %>%
        filter(permutation == x) -> haploblocks
      
      registerDoParallel(cores=nthreads)
      window_inversion<-foreach(i=1:nrow(haploblocks), .combine='rbind', .errorhandling='stop') %dopar% { 
        mutate(select(filter(convergent_cluster,
                             chromosome == as.character(haploblocks[i,1]) &
                               (Taxa1 == as.character(haploblocks[i,4]) | Taxa2 == as.character(haploblocks[i,4])) &
                               ((start >= as.numeric(haploblocks[i,2]) & start <= as.numeric(haploblocks[i,3])) | 
                                  (end >= as.numeric(haploblocks[i,2]) & end <= as.numeric(haploblocks[i,3])) | 
                                  (start < as.numeric(haploblocks[i,2]) & end > as.numeric(haploblocks[i,3])))),
                      chromosome,start,end,comparison),
               inversion=as.character(haploblocks[i,6]),
               inversion_start=as.numeric(haploblocks[i,2]),
               inversion_end=as.numeric(haploblocks[i,3]))}
      
      
      if (nrow(window_inversion) > 0) {
      window_inversion %>%
        mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
        mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
        mutate(overlap_size=(overlap_end-overlap_start)+1) %>%
        mutate(size=(end-start)+1) %>% 
        select(-inversion_start,-inversion_end) %>%
        mutate(cluster=paste0(chromosome,":",start,"-",end)) %>%
        mutate(selection_variable=paste0(comparison,"__",chromosome,"__",cluster,"__",size)) -> window_inversion
      
      window_inversion %>%
        distinct(selection_variable) %>%
        pull(selection_variable) -> loop_list
      
      rm(haploblocks)

      registerDoParallel(cores=nthreads)
      overlap_length_data<-foreach(i=1:length(loop_list), .combine='rbind', .errorhandling='stop') %dopar% { 
        reduce(IRanges(pull(filter(window_inversion,selection_variable==loop_list[i]),overlap_start),pull(filter(window_inversion,selection_variable==loop_list[i]),overlap_end))) %>%
          data.frame() %>%
          summarise(overlap_length=sum(width)) %>%
          mutate(selection=loop_list[i]) %>%
          separate(selection,into = c("comparison","chromosome","cluster","size"),sep = "__",remove = T) %>% 
          select(comparison,chromosome,cluster,size,overlap_length)}
      
      rm(window_inversion)
      
      convergent_cluster %>%
        mutate(cluster = paste0(chromosome,":",start,"-",end)) %>%
        filter(!cluster %in% pull(overlap_length_data,cluster)) %>%
        mutate(size = (as.numeric(end)-as.numeric(start))+1,
               overlap_length = 0) %>%
        select(comparison,
               chromosome,
               cluster,
               size,
               overlap_length) %>%
        rbind(.,
              overlap_length_data) -> overlap_length_data
      
      overlap_length_data %>%
        group_by(comparison) %>%
        summarise(N_cluster=n_distinct(cluster),
                  total_cluster_length=sum(as.numeric(size)),
                  total_overlap_length=sum(as.numeric(overlap_length))) %>%
        ungroup(.) %>% 
        full_join(.,
                  overlap_length_data %>%
                    filter(as.numeric(overlap_length) > 0) %>%
                    group_by(comparison) %>%
                    summarise(N_cluster_overlap=n_distinct(cluster)) %>%
                    ungroup(.)) %>%
        mutate(N_cluster_overlap=replace_na(N_cluster_overlap,0)) %>%
        mutate(porportion_number_cluster_overlap_inversion=((as.numeric(N_cluster_overlap)/as.numeric(N_cluster)))) %>%
        mutate(porportion_length_cluster_overlap_inversion=((as.numeric(total_overlap_length)/as.numeric(total_cluster_length)))) %>%
        mutate(permutation = x) %>%
        rbind(.,
              perm_result) -> perm_result
      
      rm(overlap_length_data)} else {
        
        data.frame(comparison = COMPARISON,
                   N_cluster = nrow(convergent_cluster),
                   total_cluster_length = convergent_cluster %>% mutate(size=(as.numeric(end)-as.numeric(start))+1) %>% group_by(comparison) %>% summarise(total_length=sum(as.numeric(size))) %>% pull(total_length),
                   total_overlap_length = 0,
                   N_cluster_overlap = 0,
                   porportion_number_cluster_overlap_inversion = 0,
                   porportion_length_cluster_overlap_inversion = 0,
                   permutation = x) %>%
          rbind(.,
                perm_result) -> perm_result
        rm(haploblocks,
           window_inversion)
        
      }
      
    }



fread(FILE,sep = "\t",header = T) -> overlap_length_data
  
  
  overlap_length_data %>%
  group_by(comparison) %>%
  summarise(N_cluster=n_distinct(CRA),
            total_cluster_length=sum(as.numeric(CRA_size)),
            total_overlap_length=sum(as.numeric(overlap_length))) %>%
  ungroup(.) %>%
  full_join(.,
            overlap_length_data %>%
              filter(as.numeric(overlap_length) > 0) %>%
              group_by(comparison) %>%
              summarise(N_cluster_overlap=n_distinct(CRA)) %>%
              ungroup(.)) %>%
  mutate(N_cluster_overlap=replace_na(N_cluster_overlap,0)) %>%
  mutate(porportion_number_cluster_overlap_inversion=((as.numeric(N_cluster_overlap)/as.numeric(N_cluster)))) %>%
  mutate(porportion_length_cluster_overlap_inversion=((as.numeric(total_overlap_length)/as.numeric(total_cluster_length)))) -> convergent_cluster_inversion_overlap_direction
  

  
  number_cluster_larger <- nrow(filter(perm_result,porportion_number_cluster_overlap_inversion > as.numeric(convergent_cluster_inversion_overlap_direction[1,6])))
  length_cluster_overlap <- nrow(filter(perm_result,porportion_length_cluster_overlap_inversion > as.numeric(convergent_cluster_inversion_overlap_direction[1,7])))
  permutation_count <- nrow(perm_result)
  
  convergent_cluster_inversion_overlap_direction$porportion_number_cluster_overlap_inversion_Pvalue <- number_cluster_larger/permutation_count
  convergent_cluster_inversion_overlap_direction$porportion_length_cluster_overlap_inversion_Pvalue <- length_cluster_overlap/permutation_count
  

  fwrite(convergent_cluster_inversion_overlap_direction,
         paste0(SAVE_DIR,
                "/",
                basename(FILE),
                "_Pvalue"),
         sep = "\t",
         col.names = T,
         quote = F)
