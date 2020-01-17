library(tidyverse)
library(Biostrings)

mag0 <- read_tsv('BR32_old_all.txt', 
                col_names = c('qseqid','sseqid','pident','length','mismatch',
                              'gapopen','qstart','qend','sstart','send','evalue',
                              'bitscore'))
### filter homologs > 90
mag0 %>% filter(pident > 80, evalue < 1e-20,  #by == 50?
               bitscore > 50, qseqid != sseqid) -> mag0



### did not filter out paralogous, great
write_csv(mag0, 'BR32_old_all.csv')

library(tidyverse)
library(Biostrings)
mag <- read_tsv('Mag_TPS_blast.txt',     
                col_names = c('qseqid','sseqid','pident','length','mismatch',
                              'gapopen','qstart','qend','sstart','send','evalue',
                              'bitscore')) 

### filter homologs  
mag %>% filter(pident > 80, evalue < 1e-20,
               bitscore > 50, qseqid != sseqid) -> mag2
#write_tsv(mag2, 'Mag_TPS_filtered.txt')
write_csv(mag2, 'Mag_TPS_filtered.csv')#### manually fix name 87-120(_FR13)_Oryza_sativa






library(tidyverse)
library(Biostrings)
mag_fixed <- read_csv('BR32_old_all.csv')
separate(
  mag_fixed, qseqid, c('qseqid1', 'qseqid2', 'qseqid3', 'qseqid4', 'qseqid5'),
  sep = '_', remove = FALSE) -> mag3
for (i in 1:nrow(mag3)){
  if (is.na(mag3$qseqid5[i])){
    mag3$qseqid5[i] <-  mag3$qseqid4[i]
    mag3$qseqid4[i] <-  mag3$qseqid3[i]
    mag3$qseqid3[i] <-  mag3$qseqid2[i]
    mag3$qseqid2[i] <-  NA
  }
}
for (i in 1:nrow(mag3)){
  if (!is.na(mag3$qseqid2[i])){
    mag3$qseqid1[i] <- paste(mag3$qseqid1[i], mag3$qseqid2[i], sep = '_')
  }
}
mag3$qseqid2 <- paste(mag3$qseqid3, mag3$qseqid4, sep = '_')
separate(
  mag3, sseqid, c('sseqid1', 'sseqid2', 'sseqid3', 'sseqid4', 'sseqid5'),
  sep = '_', remove = FALSE) -> mag4
mag4$sseqid2 <- paste(mag4$sseqid2, mag4$sseqid3, sep = '_')
mag4[mag4$qseqid2 == 'Lolium_perenn',]$qseqid2 -> 'Lolium_perenne' 
mag4$sseqid3 <- paste(mag4$qseqid1, mag4$sseqid1, sep = '_')
mag5 <- mag4[!duplicated(mag4$qseqid),]
mag5 %>% mutate(qseqid6 = paste(qseqid2, qseqid1, sep = '_')) -> mag6
mag6 %>% mutate(qseqid7 = paste(qseqid6, sseqid1, sep = '_')) -> mag7
write_csv(mag7, 'BR32_old_all_cleaned.csv')
## retain higher bitscore by delated duplicated ones
#####Mag_TPS.csv is the output after cleaning



MT <- read_csv('BR32_old_all_cleaned.csv')
### compute similarity
#t_id <- c('TPS1','TPS2','TPS3','TPS4','TPS5','TPS6','TPS7')

t_id <- c('tps1','tps2','tps3','tps4','tps5')


 MT %>% group_by(sseqid4) %>%    #this is not in plot
    summarise(Q5 = quantile(pident, 0.05),
              Q10 = quantile(pident, 0.05),
              Q25 = quantile(pident, 0.25),
              Q50 = quantile(pident, 0.50),
              Q75 = quantile(pident, 0.75)) 
 
 

for(i in 1:5){
     assign(paste("p", i, sep=""), 
     MT %>% filter(sseqid4 == t_id[i]) %>% 
     ggplot() +
     aes(x = pident) +
     geom_histogram(bins = 20) +
     labs(title = t_id[i])) 
}

library(egg)
plt <- ggarrange(p1, p2, p3, p4, p5, ncol=3,nrow=2) 
# plot 'TPSidentDistribu.png' saved

### compute TPS existence (unique [0/1])
TPS_n <- rep(NA, 5)
TPS_non1 <- as.list(rep(NA,5))
TPS_non2 <- as.list(rep(NA,5))

T_ID <- c('tps1','tps2','tps3','tps4','tps5')

m = length(T_ID)
for (i in 1:5){
  length(unique(filter(MT, sseqid4 == T_ID[i])$qseqid2)) -> TPS_n[i]
  unique(filter(MT, sseqid4 == T_ID[i])$qseqid2)
  unique(MT$qseqid2)[!unique(MT$qseqid2) %in% unique(filter(MT, sseqid4 == T_ID[i])$qseqid2)] -> TPS_non1[[i]]
  unique(MT$qseqid6)[!unique(MT$qseqid6) %in% unique(filter(MT, sseqid4 == T_ID[i])$qseqid6)] -> TPS_non2[[i]]
  }
TPS_existence <- data_frame(TPS = T_ID, 
                            `n/52` = TPS_n)

library(gridExtra)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(TPS_existence, rows=NULL, theme=tt)
grid.arrange(tbl,
             as.table=TRUE) ##plot 'TPS_existence'saved
TPS_non1



TPS_non2   ###this is useful




### STATs&Plots
MT %>% mutate(qseqid6 = paste(qseqid1, qseqid2, sep = '_')) -> MT
strain_id <- unique(MT$qseqid6)
host_id <- unique(MT$qseqid2)
host_occur <- rep(NA, 18)
strain_in_host <- as.list(rep(NA, 18))
for (i in 1:18){
  host_occur[i] <- 
    length(unique(filter(MT, qseqid2 == host_id[i])$qseqid6))
  strain_in_host[[i]] <- unique(filter(MT, qseqid2 == host_id[i])$qseqid6)
}
host_occrence <- tibble(host = host_id, occurence = host_occur)
host_occrence %>% 
  arrange(desc(occurence)) -> host_oc
tt2 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl2 <- tableGrob(host_oc, rows=NULL, theme=tt2)
grid.arrange(tbl2,
             as.table=TRUE)   #saved as 'host_group'
strain_in_host
#write_lines(strain_in_host, 'strain_in_host.txt')


###Overall mean
MT %>% 
  group_by(qseqid2) %>% 
  summarise(mean_pident = mean(pident)) %>% 
  arrange(desc(mean_pident)) -> TPS_mpident
tt0 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl0 <- tableGrob(TPS_mpident, rows=NULL, theme=tt0)
grid.arrange(tbl0,
             as.table=TRUE) #plot 'mean_pident_species.png'saved

###Overall median
MT %>% 
  group_by(qseqid2) %>% 
  summarise(median_pident = median(pident)) %>% 
  arrange(desc(median_pident)) -> TPS_mpident1
tt1 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1 <- tableGrob(TPS_mpident1, rows=NULL, theme=tt1)
grid.arrange(tbl1,
             as.table=TRUE) #plot 'median_pident_species.png'saved



###TPS1 median
MT %>% 
  filter(sseqid1 == 'TPS1') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS1 = median(pident)) %>% 
  arrange(desc(MedPidTPS1)) ->  TPS1
tt1_1 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_1 <- tableGrob(TPS1, rows=NULL, theme=tt1_1)
grid.arrange(tbl1_1,
             as.table=TRUE) #plot 'MedPidTPS1.png'saved

###TPS2 median
MT %>% 
  filter(sseqid1 == 'TPS2') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS2 = median(pident)) %>% 
  arrange(desc(MedPidTPS2)) ->  TPS2
tt1_2 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_2 <- tableGrob(TPS2, rows=NULL, theme=tt1_2)
grid.arrange(tbl1_2,
             as.table=TRUE) #plot 'MedPidTPS2.png'saved

###TPS3 median
MT %>% 
  filter(sseqid1 == 'TPS3') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS3 = median(pident)) %>% 
  arrange(desc(MedPidTPS3)) ->  TPS3
tt1_3 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_3 <- tableGrob(TPS3, rows=NULL, theme=tt1_3)
grid.arrange(tbl1_3,
             as.table=TRUE) #plot 'MedPidTPS3.png'saved

###TPS4 median
MT %>% 
  filter(sseqid1 == 'TPS4') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS4 = median(pident)) %>% 
  arrange(desc(MedPidTPS4)) ->  TPS4
tt1_4 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_4 <- tableGrob(TPS4, rows=NULL, theme=tt1_4)
grid.arrange(tbl1_4,
             as.table=TRUE) #plot 'MedPidTPS4.png'saved

###TPS5 median
MT %>% 
  filter(sseqid1 == 'TPS5') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS5 = median(pident)) %>% 
  arrange(desc(MedPidTPS5)) ->  TPS5
tt1_5 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_5 <- tableGrob(TPS5, rows=NULL, theme=tt1_5)
grid.arrange(tbl1_5,
             as.table=TRUE) #plot 'MedPidTPS5.png'saved


###TPS6 median
MT %>% 
  filter(sseqid1 == 'TPS6') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS6 = median(pident)) %>% 
  arrange(desc(MedPidTPS6)) ->  TPS6
tt1_6 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_6 <- tableGrob(TPS6, rows=NULL, theme=tt1_6)
grid.arrange(tbl1_6,
             as.table=TRUE) #plot 'MedPidTPS6.png'saved

###TPS7 median
MT %>% 
  filter(sseqid1 == 'TPS7') %>% 
  group_by(qseqid2) %>% 
  summarise(MedPidTPS7 = median(pident)) %>% 
  arrange(desc(MedPidTPS7)) ->  TPS7
tt1_7 <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl1_7 <- tableGrob(TPS7, rows=NULL, theme=tt1_7)
grid.arrange(tbl1_7,
             as.table=TRUE) #plot 'MedPidTPS5.png'saved
grid.arrange(tbl1_1, tbl1_2, tbl1_3, tbl1_4, tbl1_5, tbl1_6, tbl1_7, tbl1, tbl, ncol = 5, nrow = 2,
             as.table=TRUE) #plot 'MedianPid.png' saved




###extract_TPS
aa_cat <- readAAStringSet('cat.faa')
file_name <- c('tps1_extracted.fasta','tps2_extracted.fasta','tps3_extracted.fasta',
               'tps4_extracted.fasta', 'tps5_extracted.fasta')
for(i in 1:5){
  writeXStringSet(aa_cat[filter(MT,sseqid4 == 
      t_id[i], qseqid2 %in% c('Oryza_sativa', 'Triticum_aestivum'))
      $qseqid], file_name[i]) 
}

aa_cat["jgi|Magor1|10436|rna_MGG_14722T0"]



tps_extracted <- readAAStringSet('tps_extracted.fasta')




###Check duplicated events
MT %>%  mutate(qseqid8 =
                 paste(qseqid1, sseqid4, sep = '_'),
               qseqid9 = 
                 paste(qseqid2, qseqid8, sep = '_')) -> MT
  


for (i in 1:5) {
assign(paste('dup', i, sep = ''), filter(MT,sseqid4 == t_id[i])$qseqid9[
  duplicated(filter(MT,sseqid4 == t_id[i])$qseqid8)])
}


tps_dup <- c(dup1, dup2, dup3, dup4, dup5)
#filter(MT, sseqid3 %in% tps_dup)$qseqid








































