#damage.type <- "PP"

#damage.types <- c('PP', 'CPD')
#altDNAData.path <- paste("data/", altDNACode, "_data.csv", sep = "")
#UVData.path <- paste("data/", damage.type, "_all-damage-info.csv", sep="")
#UVData.folder = "data/"; UVData.filename = "_all-damage-info.csv"

UVDamageAtaltDNACounts <- function(altDNACode,
                               damage.types = c('PP', 'CPD'),
                               altDNAData.path = paste("data/", altDNACode, "/", altDNACode, "_data.csv", sep = ""), 
                               UVData.folder = "data/", UVData.filename = "_all-damage-info.csv",
                               saveTab = TRUE, saveTab.path = paste("data/", altDNACode, "/uv-damage-at-", altDNACode, "-counts.csv", sep="")){
  
  altDNAData <- data.table::fread(altDNAData.path)
  
  rtn <- data.frame(Chr = altDNAData$Chr,
                    Strand = altDNAData$Strand,
                    Start = altDNAData$Start,
                    End = altDNAData$End,
                    Quality = altDNAData$Quality)
  
  for (damage.type in damage.types){
    
    UVData.path = paste(UVData.folder, damage.type, UVData.filename, sep="")
    UVData <- data.table::fread(UVData.path)
    
    for (pdimer in substr(generateCat(ALL = FALSE, damage.type = damage.type), 1, 2)){
      
      damageOverlapCounts <- countDamageOverlap(damageData = UVData[which(UVData$damageType == paste(pdimer, damage.type, sep=".")), ],
                                                alternativeDNAData = altDNAData)
      rtn <- eval(parse(text = paste("cbind(rtn,",
                                     damage.type, ".", pdimer, ".sense.altDNA = altDNAData$", pdimer, ", ",
                                     damage.type, ".", pdimer, ".sense.overlap = damageOverlapCounts$sense,",
                                     damage.type, ".", pdimer, ".antisense.altDNA = altDNAData$", pdimer, ".rev, ",
                                     damage.type, ".", pdimer, ".antisense.overlap = damageOverlapCounts$antisense",
                                     ")", sep="")))
      
    }
    
  }
  
  if (saveTab == TRUE){write.csv(rtn, file = saveTab.path, row.names = FALSE)}
  return(rtn)
  
}

