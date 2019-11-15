library(Biostrings)
library(tidyverse)
dat_seq <- readDNAStringSet('726_v7.fasta')

dat_seq_core <- dat_seq[sort(c(2, 10, 12, 27, 46,
                               78, 98, 103, 121, 362,
                               18, 20, 28, 31, 42, 55,
                               60, 61, 86, 133, 11, 29, 
                               37, 56, 84, 13, 23, 26, 15,
                               22, 32, 63, 8, 14, 16, 35,
                               49, 53, 54, 360, 4, 5,
                               25, 40, 47, 50, 1, 21, 34, 
                               36, 43, 6, 17, 19, 39, 41,
                               68, 74, 458, 7, 9, 33, 38,
                               3, 24, 30, 100))]

dat_seq_unmapped <- dat_seq[-sort(c(2, 10, 12, 27, 46,
                                   78, 98, 103, 121, 362,
                                   18, 20, 28, 31, 42, 55,
                                   60, 61, 86, 133, 11, 29, 
                                   37, 56, 84, 13, 23, 26, 15,
                                   22, 32, 63, 8, 14, 16, 35,
                                   49, 53, 54, 360, 4, 5,
                                   25, 40, 47, 50, 1, 21, 34, 
                                   36, 43, 6, 17, 19, 39, 41,
                                   68, 74, 458, 7, 9, 33, 38,
                                   3, 24, 30, 100))]

names(dat_seq_core) <- paste('Core_Scaffold', 1:67, sep = '_')
names(dat_seq_unmapped) <- paste('Unmapped_Scaffold', 1:(911-67), sep = '_')

sum(width(dat_seq_core))
sum(width(dat_seq_unmapped))

writeXStringSet(dat_seq_core, "726_v8.fasta")
writeXStringSet(dat_seq_unmapped, "726_v8.fasta", append = T)

dat_seq2 <- readDNAStringSet('726_v8.fasta')
dat_seq2[67+74]



