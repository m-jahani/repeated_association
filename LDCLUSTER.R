################################################################################################################################################################
#################################################################LD base window clustering################################################################################
################################################################################################################################################################

#################################################################Softwares######################################################################################

#PLINK

##################################################################Arguments#########################################################################################
args = commandArgs(trailingOnly = TRUE)
NULL_FDR = args[1] # topcan.R scripr result. convergent windows that pass the FDR test. (topcand.R script result)
data = args[2] #association result (GWAS,BYPASS, or correlation), includes for columns (with out header) as follow: 1. species  2. Variable 3.SNP_ID (example Ha412HOChr01:180057) 4.window ID (example Ha412HOChr01:180001-185000) 5.asociation value (P,BF or r)
VCF = args[3] #The VCF (does not accept compress)
Threshold= args[4]#LD_threshold
Gmap=args[5] #Genetic map
threshold_cM = args[6] #Genetic Distance cut-off (5 cM) 
save_dir <- args[7]
nthreads = as.numeric(args[8])
#Runs as Rscript LDCLUSTER.R Annuus_Argophyllus_NFFD_nullw_result annuus_GWAS_result.txt annuus.vcf Annuus_threshold.90 annuus.geneticmap 5 ~ 150

#################################################################Packages#######################################################################################
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)

################################################################Output######################################################################################

#$argument1.clustering_result #indicates the position of each LD cluster and the number of convergent windows which formed the cluster

##############################################################Data preparation##################################################################################

fread(NULL_FDR,
      header  =  T) %>%
  select(Window) %>%
  distinct(Window) -> WIN # convergent windows

fread(data,
      header  =  F,
      showProgress = F) %>% 
  filter(V4 %in% pull(WIN,
                      Window)) %>%
  select(snp_id = V3,
         window = V4) -> CHR_POS #convergent SNPs ID and Window

WIN %>%
  separate(Window,
           into=c("chr",
                  "Start_End"),
           sep = ":",
           remove = T) %>%
  distinct(chr) %>%
  pull(chr) -> CHR_LIST # list of chromosomes that convergent windows are located in

fread(Threshold,
      header  =  F,
      showProgress = F) %>% 
  rename(P_distance = V1 ,
         LD_threshold = V2) -> THRESH # Ld thresholds base on physical distance

WIN %>%
  separate(Window,
           into  =  c("chr",
                      "start_end"),
           sep = ":",
           remove  =  F) %>%
  separate(start_end,
           into  =  c("start","end"),
           sep = "-", remove  =  T) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  mutate(pos = (((end-start)/2)+start)) %>%
  select(Window,
         chr,
         pos) -> conv_pos # position of convergent windows base on (((end-start)/2)+start)



fread(Gmap,
      header  =  T,
      showProgress = F) -> map #genetic map

map %>%
  group_by(chr) %>%
  arrange(pos,
          .by_group  =  TRUE) %>%
  mutate(delta_cM = cM-lag(cM,
                           default  =  first(cM)),
         delta_bp = pos-lag(pos,
                            default  =  first(pos)),
         start = lag(pos)) %>% 
  ungroup() %>%
  mutate(recomb_rate = delta_cM/delta_bp) %>%
  select(chr,
         start,
         end = pos,
         recomb_rate) %>%
  filter(!is.na(start)) %>%
  left_join(.,
            rename(map,
                   start = pos)) -> recombination # calculate recombination base on genetic map
rm(map)  
##############################################################LD calculation and Genetic map##################################################################################

registerDoParallel(cores=nthreads) 
#finds cM position and recombination rate for each convergent window base on overlapping markerd in genetic map
window_recombination <- foreach(row = 1:nrow(recombination),
                                .combine = 'rbind',
                                .errorhandling = 'stop') %dopar% { 
                                  mutate(select(filter(conv_pos,
                                                       chr  ==  as.character(recombination[row,1])  &
                                                         pos >=  as.numeric(recombination[row,2])&
                                                         pos <=  as.numeric(recombination[row,3])),
                                                Window,pos),
                                         start = as.numeric(recombination[row,2]),
                                         recomb_rate = as.numeric(recombination[row,4]),
                                         cM = as.numeric(recombination[row,5]))}
rm(conv_pos,recombination)

#estimate position of convergent windows between genetic map markers
window_recombination %>% 
  mutate(G_position = (((as.numeric(pos)-as.numeric(start))*as.numeric(recomb_rate))+as.numeric(cM))) %>%
  select(Window,G_position) -> window_recombination


registerDoParallel(cores=nthreads) 
foreach(CHR = 1:length(CHR_LIST),
        .combine = 'rbind',
        .errorhandling = 'stop') %dopar% {
          
          CHR_POS %>%
            filter(grepl(
              CHR_LIST[CHR],
              snp_id)) %>%
            select(snp_id) %>%
            fwrite(.,
                   file  =  paste0(save_dir,
                                   "/",
                                   basename(NULL_FDR),
                                   "_",
                                   CHR_LIST[CHR]),
                   sep  =  "\t",
                   col.names  =  F) # extract SNP IDs of convergent window on $CHR
          
          system(paste0("plink --vcf ",
                        VCF,
                        " --set-missing-var-ids @:# --extract ",
                        save_dir,
                        "/",
                        basename(NULL_FDR),
                        "_",
                        CHR_LIST[CHR],
                        " --r2 --ld-window-kb 999999999 --ld-window 999999999 --ld-window-r2 0 --allow-extra-chr --double-id --out ",
                        save_dir,
                        "/",
                        basename(NULL_FDR),
                        "_",
                        CHR_LIST[CHR]
          )
          )# calculate paieize LD between all the extracted SNPs
          
          system(paste0(
            "rm ",
            save_dir,
            "/",
            basename(NULL_FDR),
            "_",
            CHR_LIST[CHR])
          ) #delete the extracted SNP IDs
          
          
          fread(paste0(save_dir,
                       "/",
                       basename(NULL_FDR),
                       "_",
                       CHR_LIST[CHR],
                       ".ld"),
                header  =  T,
                showProgress  =  F) %>%
            select(SNP_A ,
                   SNP_B ,
                   R2 ) %>% #read calculated LD
            left_join(.,
                      rename(CHR_POS,
                             SNP_A=snp_id,
                             WIN_A=window)) %>% #find SNP_A and their overlapping convergent window
            filter(!is.na(WIN_A)) %>% #filter out non-overlapping SNP_As
            left_join(.,
                      rename(CHR_POS,
                             SNP_B=snp_id,
                             WIN_B=window)) %>% #find SNP_B and their overlapping convergent window
            filter(!is.na(WIN_B)) %>% #filter out non-overlapping SNP_Bs
            mutate(DW = paste(WIN_A,
                              WIN_B,
                              sep="_"),
                   R2 = as.numeric(R2)) %>%
            select(DW,
                   R2) %>%
            group_by(DW) %>%
            summarise(LD = mean(R2)) %>% #calculate mean LD for each convergent window pair
            ungroup() %>%
            separate(DW,
                     into  =  c("window_A",
                                "window_B"),
                     sep = "_",
                     remove  =  T) %>%
            separate(window_A,
                     into  =  c("chrom_A",
                                "start_end_A"),
                     sep = ":",
                     remove  =  F) %>%
            separate(start_end_A,
                     into  =  c("start_A","end_A"),
                     sep = "-",
                     remove  =  T) %>%
            separate(window_B,
                     into  =  c("chrom_B",
                                "start_end_B"),
                     sep = ":",
                     remove  =  F) %>%
            separate(start_end_B,
                     into  =  c("start_B",
                                "end_B"),
                     sep = "-",
                     remove  =  T) %>%
            mutate(P_distance = abs(as.numeric(start_A)-as.numeric(start_B))) %>% #calculate physical distance between two windows
            left_join(.,
                      THRESH) %>% #match threshold for physical distance
            mutate(LD_significant_ = (LD >= LD_threshold)) %>% #finds if passes LD threshold 
            left_join(.,
                      rename(window_recombination,
                             window_A = Window,
                             G_position_A = G_position)) %>% #finds genetic position for window_A
            left_join(.,
                      rename(window_recombination,
                             window_B = Window,
                             G_position_B = G_position)) %>% #finds genetic position for window_B
            mutate(G_distance = abs(as.numeric(G_position_A)-as.numeric(G_position_B))) %>% #calculates genetic distance between two window
            mutate(GD_significant = (G_distance >= threshold_cM)) %>% #finds if passes genetic distance threshold 
            mutate(LD_GD_sig = ifelse(LD_significant_ == "TRUE" & GD_significant == "TRUE",
                                      "TRUE",
                                      "FALSE")) %>% #find if passes both thresholds
            filter(P_distance > 0) %>% #filter out pairs windows with 0 distance (pir of a window with it self)
            select(window_A,
                   window_B,
                   LD_GD_sig) %>%
            fwrite(paste0(save_dir,
                          "/",
                          basename(NULL_FDR),
                          "_",
                          CHR_LIST[CHR],
                          ".LD_WINDOW"),
                   sep="\t",
                   col.names =T)
          
          
          system(paste0(
            "rm ",
            save_dir,
            "/",
            basename(NULL_FDR),
            "_",
            CHR_LIST[CHR],
            ".ld"))
          system(paste0(
            "rm ",
            save_dir,
            "/",
            basename(NULL_FDR),
            "_",
            CHR_LIST[CHR],
            ".log"))
          system(paste0(
            "rm ",
            save_dir,
            "/",
            basename(NULL_FDR),
            "_",
            CHR_LIST[CHR],
            ".nosex"))
          
          rbind(fread(paste0(save_dir,
                             "/",
                             basename(NULL_FDR),
                             "_",
                             CHR_LIST[CHR],
                             ".LD_WINDOW"),
                      header  =  T,
                      showProgress  =  F),
                select(fread(paste0(save_dir,
                                    "/",
                                    basename(NULL_FDR),
                                    "_",
                                    CHR_LIST[CHR],
                                    ".LD_WINDOW"),
                             header  =  T,
                             showProgress  =  F),
                       window_A = window_B,
                       window_B = window_A,
                       LD_GD_sig)) %>% 
            fwrite(paste0(save_dir,
                          "/",
                          basename(NULL_FDR),
                          "_",
                          CHR_LIST[CHR],
                          ".LD_WINDOW_TABLE"),
                   sep="\t",
                   col.names =T) #generate two way table, which means winA_winB is also avilable as WinB_winA
          
          system(paste0(
            "rm ",
            save_dir,
            "/",
            basename(NULL_FDR),
            "_",
            CHR_LIST[CHR],
            ".LD_WINDOW"))
        }

rm(CHR_POS,THRESH,window_recombination)

##############################################################LD clusteing##################################################################################


fread(data,
      header  =  F,
      showProgress = F) %>% 
  filter(V4 %in% pull(WIN,
                      Window)) %>%
  select(window = V4) %>%
  separate(window,
           into  =  c("chrom",
                      "start_end"),
           sep = ":",
           remove  =  F) %>%
  separate(start_end,
           into  =  c("start","end"),
           sep = "-",
           remove  =  T) %>%
  mutate(start = as.numeric(start)) %>%
  select(chrom,
         start,
         window) -> convergent_window # convergent windows with their chromosome and start location

for (CHR in 1:length(CHR_LIST)) {
  #filter for chromosome i
  convergent_window %>%
    filter(chrom == CHR_LIST[CHR]) %>%
    arrange(start) %>%
    select(-chrom,
           -start) %>%
    distinct(window)-> temporary #window IDs of convergent window on each $CHR
  
  fread(paste0(NULL_FDR,
               "_",
               CHR_LIST[CHR],
               ".LD_WINDOW_TABLE"),
        header  =  T,
        showProgress  =  F) -> con_sig #two way table of convergent windows pair(s)
  #empty cluster variable
  cluster <- NULL
  #fills the cluster variable up with LD clustering steps(a window pair per step)
  if (1<nrow(temporary)){
    for (o in 1:(nrow(temporary)-1)) {
      rbind(
        data.frame(
          window_A = as.character(temporary[o,1]),
          window_B = as.character(temporary[o+1,1]),
          cycle = o),
            cluster)-> cluster

    }#rows
    #calculates if each of the two windows is linked
    left_join(cluster,
              con_sig) %>%
      mutate(LD_GD_sig = replace_na(as.character(LD_GD_sig),
                                    "FALSE")) %>%
      arrange(cycle) %>%
      mutate(clustering = 0) -> cluster
  
  #It forms clusters
  for (p in 1:nrow(cluster)) {
    cluster[p,5] <- ifelse(cluster[p,3] == 1,1,
                           cluster[p,5] <- ifelse(cluster[p,5] ==  0 & cluster[p-1,4] == "TRUE" & cluster[p,4] == "TRUE" ,
                                                  cluster[p-1,5],
                                                  as.numeric(cluster[p-1,5])+1))
  }#cluster rows
  #save the clustering information
  cluster %>% filter(LD_GD_sig == "TRUE") %>%
    gather(windows,
           window_name,
           1:2) %>% 
    select(window_name,
           clustering) %>% distinct(window_name,
                                    .keep_all  =  T) -> true_cluster
  temporary %>%
    filter(!window %in% pull(distinct(true_cluster,
                                      window_name),
                             window_name)) %>%
    mutate(clustering = window) %>%
    select(window_name = window,
           clustering) %>%
    rbind(.,true_cluster) %>%
    mutate(cluster = group_indices(.,clustering)) %>%
    mutate(chrom = CHR_LIST[CHR]) %>%
    select(chrom,
           window_name,
           cluster) %>%
    fwrite(paste0(save_dir,
                  "/",
                  basename(NULL_FDR),
                  ".cluster"),
           sep  =  "\t",
           col.names =  F,
           append  =  T)

  
  }else{rbind(data.frame(chromosome = as.character(CHR_LIST[CHR]),
                         window = as.character(temporary[1,1]),
                         cluster = 1),cluster) %>%
      fwrite(paste0(save_dir,
                    "/",
                    basename(NULL_FDR),
                    ".cluster"),
             sep  =  "\t",
             col.names =  F,
             append  =  T)
  }
  rm(temporary)
  system(paste0(
    "rm ",
    save_dir,
    "/",
    basename(NULL_FDR),
    "_",
    CHR_LIST[CHR],
    ".LD_WINDOW_TABLE"))
}#chrom


#summarize cluster data
fread(paste0(NULL_FDR,
             ".cluster"),
      header = F) %>%
  select(chrom = V1,
         window_name = V2,
         cluster = V3) %>%
  separate(window_name,
           into  =  c("chr",
                      "start_end"),
           sep = ":",
           remove  =  F) %>%
  select(-chr) %>%
  separate(start_end,
           into  =  c("start",
                      "end"),
           sep = "-",
           remove  =  T) %>% 
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  group_by(chrom,
           cluster) %>%
  summarise(N_windows = n(),
            cluster_start = min(start),
            cluster_end = max(end)) %>%
  ungroup() %>%
  mutate(range = paste0(as.character(cluster_start),
                        ":",
                        as.character(cluster_end))) %>%
  mutate(size = (as.numeric(cluster_end)-as.numeric(cluster_start))+1) %>%
  select(chromosome = chrom,
         range,
         size,
         N_windows) %>%
  fwrite(paste0(save_dir,
                "/",
                str_replace(basename(NULL_FDR),"_nullw_result",""),
                ".clustering_result"),
         sep  =  "\t",
         col.names =  T)

system(paste0(
  "rm ",
  save_dir,
  "/",
  basename(NULL_FDR),
  ".cluster"))


print("DONE")
