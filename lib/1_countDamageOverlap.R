# damage.type <- "PP"; cats <- generateCat(ALL = FALSE, damage.type = damage.type)
# damageData <- UVData
# alternativeDNAData <- g4Data
# strands <- c('sense', 'antisense')

countDamageOverlap <- function(damageData, alternativeDNAData,
                               strands = c('sense', 'antisense')){

  ## this function is meant to take formation of alternative DNA structures
  ## and count the number of UV damage sites in the structure

  # damageData must have column for Start and End and chromosome
  # g4 data must have start and End coordinate, strand and chromosome

  rtn <- list()

  alternativeDNAData.strand <- factor(alternativeDNAData$Strand)

  for (strand in strands){

    if (strand == 'antisense'){
      levels(alternativeDNAData.strand) <-
        rev(levels(alternativeDNAData.strand))
    }

    # find the overlapping rows
    # query column: entries (rows) of UV damages found at a alternative
    #               DNA site (subject column)
    overlapIndices <-
      WhichOverlap(start.query = damageData$Start, end.query = damageData$End,
                   space.query = paste0(damageData$Chr, ".", damageData$Strand),
                   start.subject = alternativeDNAData$Start,
                   end.subject = alternativeDNAData$End,
                   space.subject = paste0(alternativeDNAData$Chr, ".",
                                         alternativeDNAData.strand))

    rtn[[eval(parse(text = paste("'", strand, "'", sep="")))]] <-
      unlist(lapply(1:dim(alternativeDNAData)[1], function(x, vec){
        length(which(vec == x))
        }, vec = overlapIndices[, 2]))
  }

  return(rtn)

}




