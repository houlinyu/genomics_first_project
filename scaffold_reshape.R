library(Biostrings)
library(tidyverse)
dat.seq <- readDNAStringSet('726_v3.fasta')
dat.seq[5]
#node2 dat.seq[2]
#node23 dat.seq[9]
#node30 dat.seq[14]
#scaffold 1,8 dat.seq[847] dat.seq[854]

dat.seq[919] <-subseq(dat.seq[2], 766201, width(dat.seq[2]))
dat.seq[920] <-subseq(dat.seq[9], 10407, width(dat.seq[9]))
dat.seq[921] <-subseq(dat.seq[14], 365757, width(dat.seq[14]))
dat.seq[922] <-subseq(dat.seq[847], 2868075, width(dat.seq[847]))
dat.seq[923] <-subseq(dat.seq[854], 560187, width(dat.seq[854]))

dat.seq[2] <- subseq(dat.seq[2], 1, 766200) 
dat.seq[9] <- subseq(dat.seq[9], 1, 10406)
dat.seq[14] <- subseq(dat.seq[14], 1, 365756)
dat.seq[847] <- subseq(dat.seq[847], 1, 2868074)
dat.seq[854] <- subseq(dat.seq[854], 1, 560186)

#dat.scaffold.filtered <- dat.scaffold[width(dat.scaffold[])>500]

dat.seq[5]

names(dat.seq) <- paste('Scaffold', 1:923, sep = '_')

dat.seq

writeXStringSet(dat.seq, "726_v4.fasta")

#writeXStringSet(dat.scaffold.filtered, "Desktop/726_v2.fasta", append = T)




