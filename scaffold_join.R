library(Biostrings)
library(tidyverse)
dat.seq <- readDNAStringSet('726_v4.fasta')

dat.seq[924] <- str_c(dat.seq[855],
                      dat.seq[915],
                      reverseComplement(dat.seq[860]))

dat.seq[925] <- str_c(dat.seq[859],
                      dat.seq[916])

dat.seq[926] <- str_c(dat.seq[37],
                      dat.seq[854],
                      dat.seq[850])

dat.seq[927] <- str_c(dat.seq[853], 
                      dat.seq[18],
                      reverseComplement(dat.seq[863]),
                      dat.seq[922])

dat.seq[928] <- str_c(dat.seq[8],
                      reverseComplement(dat.seq[9]),
                      dat.seq[919])

dat.seq[929] <- str_c(dat.seq[42],
                      dat.seq[923])

dat.seq[930] <- str_c(dat.seq[11],
                      reverseComplement(dat.seq[852]))
                           
dat.seq[931] <- str_c(dat.seq[19],
                      dat.seq[920])


dat.seq[932] <- subseq(dat.seq[47], 1, 62479)
dat.seq[933] <- subseq(dat.seq[47], 62632, width(dat.seq[47]))


dat.seq[934] <- str_c(dat.seq[932], 
                      dat.seq[17])

dat.seq[935] <- str_c(reverseComplement(dat.seq[933]), 
                      dat.seq[862])

dat_1 <- dat.seq[-c(855, 915, 860, 859, 916, 37, 854, 850, 853, 
                    18, 863, 922, 8, 9, 919, 42, 923, 11, 852, 19, 
                    920, 932, 933, 17, 862)]


names(dat_1) <- paste('Scaffold', 1:910, sep = '_')

dat_1

writeXStringSet(dat_1, "726_v5.fasta")

#writeXStringSet(dat.scaffold.filtered, "Desktop/726_v2.fasta", append = T)




