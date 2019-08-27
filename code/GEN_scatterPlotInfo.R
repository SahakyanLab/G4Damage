#uv.at.altDNA <- data.table::fread("data/g4/uv-damage-at-g4-counts.csv")

scatterPlotInfo <- function(altDNACode,
                            uv.at.altDNA,
                            damage.type, pdimer, strand,
                            bins,
                            raw.counts = FALSE){
  
  #raw counts will be the 
  #if false, will compute a percent change 
  
  #bin the types of g4 quadruplexes in terms of stability
  categories <- rep(0, dim(uv.at.altDNA)[1])
  for (i in 1:length(bins)){
    categories[which(uv.at.altDNA$Quality > bins[i])] <- i
  }
  uv.at.altDNA <- cbind(uv.at.altDNA, Category=categories)
  
  #determine the proportion of potential damage sites were actually observed to be damaged
  damage.sites <- pot.sites <- integer(length(bins))
  for (i in 1:(length(bins) + 1)){
    uv.at.altDNA.bin <- uv.at.altDNA[which(uv.at.altDNA$Category==i-1),]
    damage.sites[i] <- eval(parse(text = paste("sum(uv.at.altDNA.bin$",
                                               damage.type, ".", pdimer, ".", strand, ".overlap",
                                               ")", sep="")))
    pot.sites[i] <- eval(parse(text = paste("sum(uv.at.altDNA.bin$",
                                            damage.type, ".", pdimer, ".", strand, ".altDNA",
                                            ")", sep="")))
  }
  
  toPlot <- damage.sites/pot.sites
  
  #normalize 
  if (raw.counts == FALSE){toPlot <- (toPlot - toPlot[1])/toPlot[1]*100}
  
  rtn <- toPlot
  
}