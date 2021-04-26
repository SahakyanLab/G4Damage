# damage.type <- "PP"

# damage.types <- c('PP', 'CPD')
# altDNAData.path <- paste("data/", altDNACode, "_data.csv", sep = "")
# UVData.path <- paste("data/", damage.type, "_all-damage-info.csv", sep="")
# UVData.folder = "data/"; UVData.filename = "_all-damage-info.csv"

damageAtaltDNACounts <- function(damage.types,
                                 altDNAData.path,
                                 saveTab = TRUE,
                                 saveTab.path){

  # To count DNA damage at DNA loci e.g. G4
  # Dependencies: data.table, countDamageOverlap, WhichOverlap (TrantoRext)

  altDNAData <- fread(altDNAData.path)

  rtn <- data.frame(Chr = altDNAData$Chr,
                    Strand = altDNAData$Strand,
                    Start = altDNAData$Start,
                    End = altDNAData$End,
                    Quality = altDNAData$Quality)

  for (damage.type in damage.types){

    dmgData.path = paste0('raw_data/', damage.type, '_all-damage-info.csv')
    dmgData <- fread(dmgData.path)

    for (pdimer in
         substr(generateCat(damage.type = damage.type), 1, 2)){

      damageOverlapCounts <-
        countDamageOverlap(damageData = dmgData[which(dmgData$damageType ==
                                                        paste(pdimer,
                                                              damage.type,
                                                              sep=".")), ],
                           alternativeDNAData = altDNAData)
      rtn <-
        eval(parse(text =
                     paste0("cbind(rtn,",
                            damage.type, ".", pdimer,
                            ".sense.altDNA = altDNAData$", pdimer, ", ",
                            damage.type, ".", pdimer,
                            ".sense.overlap = damageOverlapCounts$sense,",
                            damage.type, ".", pdimer,
                            ".antisense.altDNA = altDNAData$", pdimer, ".rev, ",
                            damage.type, ".", pdimer,
                            ".antisense.overlap = damageOverlapCounts$antisense",
                            ")")))
    }

  }

  if (saveTab == TRUE){
    write.csv(rtn, file = saveTab.path,
              row.names = FALSE,
              quote = FALSE)
  }

  return(rtn)
}

