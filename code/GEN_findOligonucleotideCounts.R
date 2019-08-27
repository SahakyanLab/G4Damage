findOligonuceleotideCounts <- function(sequences,
                                      oligonucleotides = c("CC", "CT", "TC", "TT")){
  
  #sequences is a list of sequences that we want to find the counts of
  #finds the number of occurences of the oligonucleotides within each entry of the sequences
  #sequences can be a list or vector
  
  all.counts <- mapply(oligonucleotideFrequency, DNAStringSet(sequences),
                       MoreArgs = list(width=unique(nchar(oligonucleotides)), step=1,
                                       as.array=FALSE, as.prob=FALSE,
                                       fast.moving.side="right", with.labels=TRUE))
  
  pdimer.counts <- all.counts[which(row.names(all.counts) %in% oligonucleotides), ]
  pdimer.counts.rev <- all.counts[which(row.names(all.counts) %in% as.character(reverseComplement(DNAStringSet(oligonucleotides)))), ]
  
  row.names(pdimer.counts) <- oligonucleotides[order(oligonucleotides)]
  row.names(pdimer.counts.rev) <- paste(oligonucleotides[order(as.character(reverseComplement(DNAStringSet(oligonucleotides))))], "rev", sep=".")
  
  rtn <- list(counts = t(pdimer.counts),
              rev.counts = t(pdimer.counts.rev))
  
  return(rtn)
  
}