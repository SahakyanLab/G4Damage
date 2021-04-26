###############################################################################
# This function takes a vector, each element of which contains a single      ##
## sequence, and converts it into a fasta format with the header names being ##
## either the ones provided, or just the position identifiers (line numbers) ##
## (default)in the original raw sequence file.                               ##
###############################################################################
# REQUIRES make.fasta.R!!!
###############################################################################
make.multifasta <- function(RAWSEQ, UPPER=TRUE, ncharPline=60, HEADER=NULL){
  MULTIFASTA <- NULL
  for(i in 1:length(RAWSEQ)){
  
    if(is.null(HEADER[1])){
      header <- i
    } else {
      header <- HEADER[i]
    }
  
    MULTIFASTA <- c(MULTIFASTA, make.fasta(SEQ=unlist(strsplit(RAWSEQ[i],"")), 
                                 UPPER=UPPER, HEADER=header, ncharPline=ncharPline))
  }
  return(MULTIFASTA)
}
###############################################################################
