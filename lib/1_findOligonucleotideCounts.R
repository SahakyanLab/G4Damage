findOligonuceleotideCounts <- function(sequences, oligonucleotides){
  
  # Find the number of occurrences of the oligonucleotides within each entry of
  # the sequences.
  # Sequences is a <list> or <vector> of sequences that we want to find the
  # counts of.

  all.counts <- mapply(oligonucleotideFrequency, DNAStringSet(sequences),
                       MoreArgs = list(width = unique(nchar(oligonucleotides)),
                                       step = 1, as.array = FALSE,
                                       as.prob = FALSE,
                                       fast.moving.side = "right",
                                       with.labels = TRUE))
  
  pdimer.counts <- all.counts[which(row.names(all.counts) %in%
                                      oligonucleotides), ]
  pdimer.counts.rev <- all.counts[which(row.names(all.counts) %in%
                                          as.character(reverseComplement(
                                            DNAStringSet(oligonucleotides)))), ]
  
  # If damage pattern is only one, convert vector to matrix
  if (is.null(dim(pdimer.counts))) {
    pdimer.counts <- matrix(pdimer.counts, nrow = 1)
    pdimer.counts.rev <- matrix(pdimer.counts.rev, nrow = 1)
  }

  row.names(pdimer.counts) <- oligonucleotides[order(oligonucleotides)]
  row.names(pdimer.counts.rev) <- paste(oligonucleotides[order(
    as.character(reverseComplement(DNAStringSet(oligonucleotides))))],
    "rev", sep = ".")
  
  rtn <- list(counts = t(pdimer.counts),
              rev.counts = t(pdimer.counts.rev))
  
  return(rtn)
  
}