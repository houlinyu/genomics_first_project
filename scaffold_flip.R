library(Biostrings)
library(tidyverse)
dat_seq <- readDNAStringSet('726_v5.fasta')



dat_seq[911] <- subseq(dat_seq[5], 1, 387678)
dat_seq[912] <- subseq(dat_seq[5], 387679, 474951)


dat_seq[914] <- str_c(dat_seq[911],
                      reverseComplement(dat_seq[912]),
                      dat_seq[913])

dat_1 <- dat_seq[-c(5, 911, 912, 913)]
dat_2 <- dat_1[order(desc(width(dat_1)))]

names(dat_2) <- paste('Scaffold', 1:910, sep = '_')

writeXStringSet(dat_2, "726_v6.fasta")


#writeXStringSet(dat.scaffold.filtered, "Desktop/726_v2.fasta", append = T)




