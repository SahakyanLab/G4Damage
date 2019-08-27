## FUNCTION ####################################################################
GenomeGetKmers <- function(k = 11,
                          nproc = 1,
                          PATH.genome = "human_genome_unmasked_37.73/", 
                          # "/Users/sirius/Database/human_genome_unmasked_37.73/"
                          genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
                          genome.numer  = c(1:22,"X","Y"),
                          fastafile.ending  = ".fa"){

  registerDoParallel( makeCluster(nproc) )
  ##------------------------
  foreach(i = genome.numer) %dopar% {
    print(i, quote=FALSE)
    fasta <- readfasta(filename = paste(PATH.genome,"/",genome.prefix,i,fastafile.ending,sep=""),
                       split = TRUE)
    kmers.count <- getKmers(seq.string=fasta$seq, k=k, method="Biostrings")#changed from method="index"
    write.table(as.table(kmers.count), 
                file=paste("chr_",i,"_kmer_",k,"_tbl.txt",sep=""),
                row.names=FALSE)
  }
  ##------------------------

}
## FUNCTION ####################################################################

# EXAMPLE:
#
#TrantoR.lib = "TrantoR/" # "/Users/sirius/GIT/TrantoR/"
#source(paste(TrantoR.lib,"/GEN_readfasta.R", sep=""))
#source(paste(TrantoR.lib,"/GEN_Seq2Index.R", sep=""))
#source(paste(TrantoR.lib,"/GEN_getKmers.R", sep=""))
#source(paste(TrantoR.lib,"/GEN_GenomeGetKmers.R", sep=""))
#
#library(compiler)
#library(foreach)
#library(gtools)
#library(doMC)
#
#readfasta      <- cmpfun(readfasta)
#Seq2Index      <- cmpfun(Seq2Index)
#getKmers       <- cmpfun(getKmers)
#GenomeGetKmers <- cmpfun(GenomeGetKmers)
#
#GenomeGetKmers(k = 7,
#              nproc = 1,
#              PATH.genome = "human_genome_unmasked_37.73/", 
#              # "/Users/sirius/Database/human_genome_unmasked_37.73/"
#              genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
#              genome.numer  = c(1:22,"X","Y"),
#              fastafile.ending  = ".fa")
## Human chr 1 takes ~ 5.4 GB RAM
#


