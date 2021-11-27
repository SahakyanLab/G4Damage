# damage.type <- "PP"

# damage.types <- c('PP', 'CPD')
# altDNAData.path <- paste("data/", altDNACode, "_data.csv", sep = "")
# UVData.path <- paste("data/", damage.type, "_all-damage-info.csv", sep="")
# UVData.folder = "data/"; UVData.filename = "_all-damage-info.csv"

damageAtaltDNACounts <- function(damage.types,
                                 altDNAData.path,
                                 strand.sensitive,
                                 include.weight,
                                 saveTab.path = NULL){
  
  # To count DNA damage at DNA loci e.g. G4
  # Dependencies: data.table, countDamageOverlap, WhichOverlap (TrantoRext)
  
  chrs <- fread(altDNAData.path, showProgress = FALSE)[, unique(Chr)]
  
  rtn <- lapply(chrs, function(chr) {
    print(chr)
    
    altDNAData <- fread(altDNAData.path, showProgress = FALSE, nThread = 1)[Chr == chr]
    if (!strand.sensitive) altDNAData[, Strand := "+"]
    
    rtn <- altDNAData[, .(Chr, Strand, Start, End, Quality)]
    
    for (damage.type in damage.types){
      print(damage.type)
      
      dmgData.path = paste0('raw_data/', damage.type, "/", chr, ".csv.gz")
      if (!file.exists(dmgData.path)) next
      
      dmgData <- fread(dmgData.path, showProgress = FALSE, nThread = 1)
      if (nrow(dmgData) == 0) next
      
      for (pdimer in
           substr(generateCat(damage.type = damage.type), 1, 2)){
        
        # For single base damage e.g. 8-oxoG
        pdimer <- sub("\\.$", "", pdimer)
        
        damageOverlapCounts <- countDamageOverlap(
          damageData = dmgData[damageType == paste0(pdimer, ".", damage.type), ],
          alternativeDNAData = altDNAData,
          strands = if(!strand.sensitive) "sense" else c('sense', 'antisense'),
          include.weight = include.weight)

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
    rtn
  })
  rtn <- rbindlist(rtn, fill = TRUE)
  
  if (!is.null(saveTab.path)) fwrite(rtn, saveTab.path)
  
  return(rtn)
}

