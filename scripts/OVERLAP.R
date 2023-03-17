library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
library(IRanges)

options(sipen = 9999)


##################################################################Arguments#########################################################################################
args = commandArgs(trailingOnly = TRUE)
HAPLOBLOCKS = args[1] 
CRAS = args[2] 
COMPARISON = args[3] 
SAVE_DIR= args[4]
nthreads = as.numeric(args[5])
#####################################################################################################################################################################


fread(HAPLOBLOCKS,
      header = T) -> haploblocks.bed

read.table(CRAS,
           sep = "\t",
           header = T) %>%
  mutate(direction = COMPARISON) %>%
  filter(chromosome != "chromosome") %>%
  distinct(chromosome,range,size,N_windows,direction) %>%
  separate(direction, into = c("Taxa1","Taxa2"),sep = "_",remove = F) %>%
  separate(range,into = c("start","end"),sep = ":",remove = F) %>%
  mutate(start=as.numeric(start),end=as.numeric(end))-> convergent_cluster

#overlap between inversion and conversion regions(including partial overlap)
registerDoParallel(cores=100)
window_inversion<-foreach(i=1:nrow(haploblocks.bed), .combine='rbind', .errorhandling='stop') %dopar% { 
  mutate(select(filter(convergent_cluster,
                       chromosome == as.character(haploblocks.bed[i,1]) &
                         (Taxa1 == as.character(haploblocks.bed[i,4]) | Taxa2 == as.character(haploblocks.bed[i,4])) &
                         ((start >= as.numeric(haploblocks.bed[i,2]) & start <= as.numeric(haploblocks.bed[i,3])) | 
                            (end >= as.numeric(haploblocks.bed[i,2]) & end <= as.numeric(haploblocks.bed[i,3])) | 
                            (start < as.numeric(haploblocks.bed[i,2]) & end > as.numeric(haploblocks.bed[i,3])))),
                direction,chromosome,range,start,end,size,N_windows),
         inversion=as.character(haploblocks.bed[i,5]),
         inversion_start=as.numeric(haploblocks.bed[i,2]),
         inversion_end=as.numeric(haploblocks.bed[i,3]))}

window_inversion %>%
  mutate(overlap_start=ifelse(start>=inversion_start,start,inversion_start)) %>%
  mutate(overlap_end=ifelse(end>=inversion_end,inversion_end,end)) %>%
  mutate(overlap_size=(overlap_end-overlap_start)+1) %>%
  select(-inversion_start,-inversion_end) %>%
  mutate(cluster=paste0(chromosome,":",range)) %>%
  mutate(selection_variable=paste0(direction,"__",chromosome,"__",cluster,"__",size)) -> window_inversion


window_inversion %>%
  distinct(selection_variable) %>%
  pull(selection_variable) -> loop_list

registerDoParallel(cores=100)
overlap_length_data<-foreach(i=1:length(loop_list), .combine='rbind', .errorhandling='stop') %dopar% { 
  reduce(IRanges(pull(filter(window_inversion,selection_variable==loop_list[i]),overlap_start),pull(filter(window_inversion,selection_variable==loop_list[i]),overlap_end))) %>%
    data.frame() %>%
    summarise(overlap_length=sum(width)) %>%
    mutate(selection=loop_list[i]) %>%
    separate(selection,into = c("direction","chromosome","cluster","size"),sep = "__",remove = T) %>% 
    select(direction,chromosome,cluster,size,overlap_length)}

convergent_cluster %>%
  mutate(cluster=paste0(chromosome,":",range)) %>%
  mutate(selection_variable=paste0(direction,"__",chromosome,"__",cluster,"__",size)) %>%
  filter(!selection_variable %in% loop_list) %>%
  mutate(overlap_length=0) %>%
  select(direction,chromosome,cluster,size,overlap_length) %>%
  rbind(.,overlap_length_data) %>%
  full_join(.,
            select(window_inversion,
                   direction,
                   chromosome,
                   inversion)) %>%
  distinct(direction,
           chromosome,
           cluster,
           size,
           overlap_length,
           inversion) %>%
  group_by(direction,
           chromosome,
           cluster,
           size) %>%
  mutate(haploblock_num=paste0("inv_",row_number())) %>%
  ungroup() %>%
  spread(haploblock_num,inversion) %>%
  unite(ovelaps, starts_with("inv_"), na.rm = TRUE, remove = TRUE, sep = ",") %>%
  select(comparison = direction,
         CRA = cluster,
         CRA_size = size,
         ovelap_haaploblock = ovelaps,
         overlap_length) %>%
  mutate(ovelap_haaploblock=ifelse(overlap_length==0,NA,ovelap_haaploblock)) %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                basename(CRAS),
                "_overalp_haploblocks"),
         sep = "\t",
         col.names = T,
         quote = F,
         na = "NA")
