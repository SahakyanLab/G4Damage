#uv.at.altDNA <- data.table::fread("data/g4/uv-damage-at-g4-counts.csv")

scatterPlotInfo <- function(dmg.at.G4,
                            dmg.name, pdimer, strand,
                            bins,
                            raw.counts = FALSE){

  # raw counts will be the default
  # if false, will compute a percent change

  # bin the types of g4 quadruplexes in terms of stability
  categories <- rep(0, dim(dmg.at.G4)[1])
  for (i in 1:length(bins)){
    categories[which(dmg.at.G4$Quality > bins[i])] <- i
  }
  dmg.at.G4 <- cbind(dmg.at.G4, Category=categories)

  # Determine the proportion of potential damage sites that were actually
  # observed to be damaged
  damage.sites <- pot.sites <- integer(length(bins))
  for (i in 1:(length(bins) + 1)){
    dmg.at.G4.bin <- dmg.at.G4[which(dmg.at.G4$Category==i-1),]
    damage.sites[i] <-
      eval(parse(text = paste0("sum(dmg.at.G4.bin$", dmg.name, ".",
                               pdimer, ".", strand, ".overlap, na.rm = T)")))
    pot.sites[i] <-
      eval(parse(text = paste0("sum(dmg.at.G4.bin$", dmg.name, ".",
                               pdimer, ".", strand, ".altDNA, na.rm = T)")))
  }

  toPlot <- damage.sites/pot.sites

  #normalize
  if (raw.counts == FALSE){toPlot <- (toPlot - toPlot[1])/toPlot[1]*100}

  rtn <- toPlot

}
