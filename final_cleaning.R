####Read.fasta
library(Biostrings)
dat_seq <- readDNAStringSet('UMASS_Fom_PHW726_1.fasta')

library(tidyverse)

#writeXStringSet(dat_seq[c("Unmapped_Scaffold_250", 
#                        "Unmapped_Scaffold_327",  
#                        "Unmapped_Scaffold_441",
#                        "Unmapped_Scaffold_472")],
#                        "fo_blast.fasta")

###do not need to run
#seq_1 <- as_tibble(read.table("726_v10_pacbio_coverageperbase.txt", sep = "", header = T))
#seq_1 %>% group_by(Core_Scaffold_1) %>% summarise(S = sum(X8)) %>% filter(S == 0) -> seq_pac
#pacbionoreads_id <- seq_pac$Core_Scaffold_1
#dat_seq1 <- dat_seq[pacbionoreads_id]
#dat_seq1[width(dat_seq1) >= 1000] -> dat_seq1_1

####get_filtering_list_from_illumina
illumina <- as_tibble(read.table("726_v10_illumina_coverageperbase.txt", sep = "", header = T))
illumina %>% group_by(Core_Scaffold_1) %>% summarise(M = mean(X27)) %>% 
                         filter(M <= 11) -> dat_seq_filter_out
illulowreads_id <- dat_seq_filter_out$Core_Scaffold_1

####filter_out_below_11
dat_seq2 <- dat_seq[illulowreads_id]
dat_seq3 <- dat_seq[!dat_seq %in% dat_seq2]
dat_seq_largerthan1k <- dat_seq3[width(dat_seq3) >= 1000]
dat_seq_largerthan1k[c("Unmapped_Scaffold_365", "Unmapped_Scaffold_419", 
                       "Unmapped_Scaffold_519", "Core_Scaffold_53")] -> dat_seq4
dat_seq5 <- dat_seq_largerthan1k[!dat_seq_largerthan1k %in% dat_seq4]

####remove_adaptors
##info
#Sequence name, length, span(s), apparent source
#Core_Scaffold_13	1034886	562,252..562,278	adaptor:NGB00866.1
#Core_Scaffold_22	716004	48,083..48,109	adaptor:NGB00813.1
#Unmapped_Scaffold_143	26962	2444..2473	adaptor:NGB01064.1


dat_seq5["Core_Scaffold_13"] <- str_c(
                                      subseq(dat_seq5["Core_Scaffold_13"], 1, 562251),
                                      'NNNNNNNNNN',
                                      subseq(dat_seq5["Core_Scaffold_13"], 562279, 
                                             width(dat_seq5["Core_Scaffold_13"])))
dat_seq5["Core_Scaffold_22"] <- str_c(
                                      subseq(dat_seq5["Core_Scaffold_22"], 1, 48082),
                                      'NNNNNNNNNN',
                                      subseq(dat_seq5["Core_Scaffold_22"], 48110, 
                                             width(dat_seq5["Core_Scaffold_22"])))
dat_seq5["Unmapped_Scaffold_143"] <- str_c(
                                           subseq(dat_seq5["Unmapped_Scaffold_143"], 1, 2443),
                                           'NNNNNNNNNN',
                                           subseq(dat_seq5["Unmapped_Scaffold_143"], 2474, 
                                           width(dat_seq5["Unmapped_Scaffold_143"])))

####a_test
#subseq(dat_seq5["Core_Scaffold_13"], 562251, 562263)

dat_seq6 <- dat_seq5[1:52]

dat_seq7 <- dat_seq5[53:66]
names(dat_seq7) <- paste('Core_Scaffold', 53:66, sep = '_')

dat_seq8 <- dat_seq5[67:582]
names(dat_seq8) <- paste('Unmapped_Scaffold', 1:(582-67+1), sep = '_')

dat_seq9 <- dat_seq5[583]

dat_seq10 <- c(dat_seq6, dat_seq7, dat_seq8, dat_seq9)

writeXStringSet(dat_seq10, "final.fasta")


sum(width(dat_seq10[1:66]))

