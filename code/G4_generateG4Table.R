#G4.data.path = "data/PQSdata.txt"
#saveTab.path = "data/g4/g4_data.csv"

## ---- generateG4Table-CODE
# description

generateG4Table <- function(G4.data.path = "data/g4/PQSdata.txt", 
                            saveTab = TRUE, saveTab.path = "data/g4/g4_data.csv"){
  
  ##g4 data table should have the following columns:
  #chr: chromosome where the g-quadruplex was observed
  #genomic.start: genomic start coordinate
  #genomic.end: genomic end coordinate
  #relseq: sequence on the strand where the g-quadruplex was observed
  #mm: quality score
  #also any additional columns are fine (can be in any order)
  
  source("code/GEN_findOligonucleotideCounts.R")
  
  G4.data <- data.table::fread(G4.data.path)
  
  #get rid of quadruplexes without quality scores
  G4.data <- G4.data[which(is.na(G4.data$mm) == FALSE),]
  
  countList <- findOligonuceleotideCounts(sequences = G4.data$relseq)
  
  G4.data <- data.frame(Chr = G4.data$chr,
                        Quality = G4.data$mm,
                        Strand = G4.data$str,
                        Start = G4.data$genomic.start,
                        End = G4.data$genomic.end,
                        RelSeq = G4.data$relseq,
                        countList$counts,
                        RevComp = as.character(reverseComplement(DNAStringSet(G4.data$relseq))),
                        countList$rev.counts)
  
  if (saveTab == TRUE){write.csv(G4.data, file = saveTab.path, row.names = FALSE)}
  return(G4.data)
  
}
