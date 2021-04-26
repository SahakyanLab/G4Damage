#uv.at.altDNA <- data.table::fread("data/g4/uv-damage-at-g4-counts.csv")

scatterPlotInfo <- function(altDNACode,
                            dmg.at.altDNA,
                            damage.type, pdimer, strand,
                            bins,
                            raw.counts = FALSE){

  # raw counts will be the default
  # if false, will compute a percent change

  # bin the types of g4 quadruplexes in terms of stability
  categories <- rep(0, dim(dmg.at.altDNA)[1])
  for (i in 1:length(bins)){
    categories[which(dmg.at.altDNA$Quality > bins[i])] <- i
  }
  dmg.at.altDNA <- cbind(dmg.at.altDNA, Category=categories)

  # Determine the proportion of potential damage sites that were actually
  # observed to be damaged
  damage.sites <- pot.sites <- integer(length(bins))
  for (i in 1:(length(bins) + 1)){
    dmg.at.altDNA.bin <- dmg.at.altDNA[which(dmg.at.altDNA$Category==i-1),]
    damage.sites[i] <-
      eval(parse(text = paste0("sum(dmg.at.altDNA.bin$", damage.type, ".",
                               pdimer, ".", strand, ".overlap)")))
    pot.sites[i] <-
      eval(parse(text = paste0("sum(dmg.at.altDNA.bin$", damage.type, ".",
                               pdimer, ".", strand, ".altDNA)")))
  }

  toPlot <- damage.sites/pot.sites

  #normalize
  if (raw.counts == FALSE){toPlot <- (toPlot - toPlot[1])/toPlot[1]*100}

  rtn <- toPlot

}
