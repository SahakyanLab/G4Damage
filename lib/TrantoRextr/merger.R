genome.numer = c(1:22,"X","Y")

###############
for(k in c(2,4,6,8,10,12)){

  count <- 1
  merged.kmer <- NULL
  for(chr in genome.numer){
    print(c(k, chr))
    DB2merge <- paste("GenomeKmer_",k,"/chr_",chr,"_kmer_",k,"_tbl.txt",sep="")
    if(count == 1){
      merged.kmer <- read.table(DB2merge, header=TRUE)
    } else {
      merged.kmer$Freq <- merged.kmer$Freq + read.table(DB2merge, header=TRUE)$Freq
    }
    count <- count + 1
  }

  write.table( merged.kmer,
               row.names = FALSE,
               file=paste("GenomeKmer_",k,"/chr_ALL_kmer_",k,"_tbl.txt",sep="") )
}
###############




