generateG4Table <- function(G4.data.path, dmg.pattern,
                            saveTab.path=NULL){
  
  # AIM: Count dmg.pattern within G4 sequences and update the G4 table.
  #      Only relevant columns are retained.
  #
  # G4.data.path (G4 table) should have the following columns:
  # chr           : chromosome where the g-quadruplex was observed
  # genomic.start : genomic start coordinate
  # genomic.end   : genomic end coordinate
  # relseq        : sequence on the strand where the g-quadruplex was observed
  # mm            : quality score
  # also any additional columns are fine (can be in any order)
  #
  # Dependencies: "lib/GEN_findOligonucleotideCounts.R"
  
  G4.data <- fread(G4.data.path)

  # get rid of quadruplexes without quality scores
  G4.data <- G4.data[which(is.na(G4.data$mm) == FALSE),]
  #print(G4.data)
  
  if (!"NN" %in% dmg.pattern) {
    countList <- findOligonuceleotideCounts(sequences = G4.data$relseq,
                                            oligonucleotides = dmg.pattern)
  } else {
    countList <- list(counts = matrix(nchar(G4.data$relseq)))
    colnames(countList$counts) <- dmg.pattern
    countList$rev.counts <- countList$counts
  }
  
  G4.data <- data.table(Chr = G4.data$chr,
                        Quality = G4.data$mm,
                        Strand = G4.data$str,
                        Start = G4.data$genomic.start,
                        End = G4.data$genomic.end,
                        RelSeq = G4.data$relseq,
                        countList$counts,
                        RevComp = as.character(reverseComplement(
                          DNAStringSet(G4.data$relseq))),
                        countList$rev.counts)
  
  if (!is.null(saveTab.path)) fwrite(G4.data, saveTab.path)
  
  return(G4.data)
}
