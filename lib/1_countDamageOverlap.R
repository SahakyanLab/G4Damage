countDamageOverlap <- function(damageData, G4Data,
                               strands = c('sense', 'antisense'),
                               include.weight = FALSE){

  ## this function is meant to take formation of alternative DNA structures
  ## and count the number of UV damage sites in the structure

  # damageData must have column for Start and End and chromosome
  # g4 data must have start and End coordinate, strand and chromosome

  rtn <- list()

  G4Data.strand <- factor(G4Data$Strand)

  for (strand in strands){

    if (strand == 'antisense'){
      levels(G4Data.strand) <-
        rev(levels(G4Data.strand))
    }

    # find the overlapping rows
    # query column: entries (rows) of UV damages found at a alternative
    #               DNA site (subject column)
    overlapIndices <-
      WhichOverlap(start.query = damageData$Start, end.query = damageData$End,
                   space.query = paste0(damageData$Chr, ".", damageData$Strand),
                   start.subject = G4Data$Start,
                   end.subject = G4Data$End,
                   space.subject = paste0(G4Data$Chr, ".",
                                          G4Data.strand))#; stop()

    if (include.weight) {
      damageData <- damageData[overlapIndices[, 1]]
      damageData[, Subject := overlapIndices[, 2]]
      damageData <- damageData[, .(wCount = sum(Weight)), by = Subject]
      zeros <- damageData[, (1:nrow(G4Data))[!1:nrow(G4Data) %in% Subject]]
      damageData <- rbind(damageData, data.table(Subject = zeros, wCount = 0))
      setkey(damageData, Subject)
      rtn[[strand]] <- damageData$wCount
    } else {
      rtn[[eval(parse(text = paste("'", strand, "'", sep="")))]] <-
        unlist(lapply(1:dim(G4Data)[1], function(x, vec){
          length(which(vec == x))
        }, vec = overlapIndices[, 2]))
    }
  }

  return(rtn)

}
