TrantoR.lib = "../../TrantoR/" # "/Users/sirius/GIT/TrantoR/"
source(paste(TrantoR.lib,"/GEN_readfasta.R", sep=""))
source(paste(TrantoR.lib,"/GEN_Seq2Index.R", sep=""))
source(paste(TrantoR.lib,"/GEN_getKmers.R", sep=""))
source(paste(TrantoR.lib,"/GEN_GenomeGetKmers.R", sep=""))

library(compiler)
library(foreach)
library(gtools)
library(doParallel)

readfasta      <- cmpfun(readfasta)
Seq2Index      <- cmpfun(Seq2Index)
getKmers       <- cmpfun(getKmers)
GenomeGetKmers <- cmpfun(GenomeGetKmers)

GenomeGetKmers(k = 1,
               nproc = 5,
               PATH.genome = "../../human_genome_unmasked_37.73/", 
               # "/Users/sirius/Database/human_genome_unmasked_37.73/"
               genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
               genome.numer  = c(22:2,"X","Y",1),
               fastafile.ending  = ".fa")
# Human chr 1 takes ~ 5.2 GB RAM
              
              
