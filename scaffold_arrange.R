library(Biostrings)
library(tidyverse)
dat.seq <- readDNAStringSet('726_v6.fasta')




dat.seq[911] <-subseq(dat.seq[3], 1105601, width(dat.seq[3]))


dat.seq[3] <- subseq(dat.seq[3], 1, 1105600) 
dat_2 <- dat.seq[order(desc(width(dat.seq)))]


names(dat_2) <- paste('Scaffold', 1:911, sep = '_')
writeXStringSet(dat_2, "726_v7.fasta")




