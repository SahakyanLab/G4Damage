damageAtG4Counts <- function(dmg.names,
                             G4.path,
                             strand.sensitive,
                             include.weight,
                             saveTab.path = NULL){

  # To count DNA damage at DNA loci e.g. G4
  # Dependencies: data.table, countDamageOverlap, WhichOverlap (TrantoRext)

  chrs <- fread(G4.path, showProgress = FALSE)[, unique(Chr)]

  rtn <- lapply(chrs, function(chr) {
    print(chr)

    G4Data <- fread(G4.path, showProgress = FALSE, nThread = 1)[Chr == chr]
    if (!strand.sensitive) G4Data[, Strand := "+"]

    rtn <- G4Data[, .(Chr, Strand, Start, End, Quality)]

    for (dmg.name in dmg.names){
      print(dmg.name)

      dmgData.path = paste0('raw_data/', dmg.name, "/", chr, ".csv.gz")
      if (!file.exists(dmgData.path)) next

      dmgData <- fread(dmgData.path, showProgress = FALSE, nThread = 1)
      if (nrow(dmgData) == 0) next

      for (pdimer in
           substr(generateCat(dmg.name = dmg.name), 1, 2)){

        # For single base damage e.g. 8-oxoG
        pdimer <- sub("\\.$", "", pdimer)

        damageOverlapCounts <- countDamageOverlap(
          damageData = dmgData[damageType == paste0(pdimer, ".", dmg.name), ],
          G4Data = G4Data,
          strands = if(!strand.sensitive) "sense" else c('sense', 'antisense'),
          include.weight = include.weight)

        rtn <-
          eval(parse(text =
                       paste0("cbind(rtn,",
                              dmg.name, ".", pdimer,
                              ".sense.altDNA = G4Data$", pdimer, ", ",
                              dmg.name, ".", pdimer,
                              ".sense.overlap = damageOverlapCounts$sense,",
                              dmg.name, ".", pdimer,
                              ".antisense.altDNA = G4Data$", pdimer, ".rev, ",
                              dmg.name, ".", pdimer,
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

