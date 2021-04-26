## FUNCTION ####################################################################
Genomeregexfinder <- function(qpoutname         = "human_cPQS.qpout",
                              PATH.genome       = "/Users/sirius/Database/human_genome_unmasked_37.73/",
                              genome.prefix     = "Homo_sapiens.GRCh37.73.dna.chromosome.",
                              genome.numer      = c(1:22,"X","Y"),
                              fastafile.ending  = ".fa",
                              plus.strand.regex = "[cC]{1}([gG]{3,}[NATGCnatgc]{1,12}){3}[gG]{3,}",
                              minus.strand.regex= "([cC]{3,}[NATGCnatgc]{1,12}){3}[cC]{3,}[gG]{1}"){

  write("CHR START.POS END.POS STRAND LENGTH SEQUENCE", file=qpoutname)

  ##------------------------
  for(i in genome.numer){   
    print(i, quote=F)
    fasta <- readfasta(paste(PATH.genome,"/",genome.prefix,i,fastafile.ending,sep=""))
  
    QP.RESULTS <- DNAregexfinder(seq=fasta$seq, type="other",
                          plus.strand.regex=plus.strand.regex,
                          minus.strand.regex=minus.strand.regex, perl=TRUE)
   
    ##################                                           
    if(QP.RESULTS$num.occ>0){
   
      write( paste( rep(i, length(QP.RESULTS$start.pos)), 
                    QP.RESULTS$start.pos,
                    QP.RESULTS$start.pos + QP.RESULTS$seq.length - 1,
                    QP.RESULTS$strand,
                    QP.RESULTS$seq.length,
                    QP.RESULTS$sequence, sep=" "), file=qpoutname, append=TRUE)
    }
    ################## 
  
  }
  ##------------------------
}
## FUNCTION ####################################################################


