#uv.at.altDNA <- data.table::fread("data/g4/uv-damage-at-g4-counts.csv")

scatterPlotInfo <- function(dmg.at.G4,
                            dmg.name, pdimer, strand,
                            bins,
                            raw.counts = FALSE,
                            cal.p.vals = FALSE){

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
    dmg.at.G4.bin <- dmg.at.G4[which(dmg.at.G4$Category == i-1),]
    
    print(nrow(dmg.at.G4.bin))
    
    damage.sites[i] <-
      eval(parse(text = paste0("sum(dmg.at.G4.bin$", dmg.name, ".",
                               pdimer, ".", strand, ".overlap, na.rm = T)")))
    pot.sites[i] <-
      eval(parse(text = paste0("sum(dmg.at.G4.bin$", dmg.name, ".",
                               pdimer, ".", strand, ".altDNA, na.rm = T)")))
  }

  toPlot <- damage.sites / pot.sites

  if (cal.p.vals) {
    p.vals.1 <- sapply(seq_along(damage.sites), function(i) {
      fisher.test(matrix(c(damage.sites[1], pot.sites[1] - damage.sites[1],
                           damage.sites[i], pot.sites[i] - damage.sites[i]),
                         ncol = 2))$p.value
    })
    
    p.vals.2 <- sapply(seq_along(damage.sites), function(i) {
      zs <- (toPlot[1] - toPlot[i]) /
        sqrt(toPlot[1] * (1 - toPlot[1]) / pot.sites[1] +
               toPlot[i] * (1 - toPlot[i]) / pot.sites[i])
      pnorm(-abs(zs)) * 2
    })
    
    # [DEV] investigate p-value
    a <<- damage.sites
    b <<- pot.sites
    print(pot.sites)
    
    lr.model <- glm(cbind(damage.sites, pot.sites - damage.sites) ~
                      factor(seq_along(damage.sites)),
                    family = "binomial")
    lr.anova <- anova(lr.model, test = "Chisq")$`Pr(>Chi)`[2]
    
    cat("Fisher Test:", p.vals.1, "\n")
    cat("Z-Test:", p.vals.2, "\n")
    cat("Binomial Anova:", lr.anova, "\n")
  } else {
    p.vals.1 <- p.vals.2 <- lr.anova <- NULL
  }
  
  cnt <- signif(toPlot, 2)
  #normalize
  if (raw.counts == FALSE){toPlot <- (toPlot - toPlot[1])/toPlot[1]*100}
  
  

  rtn <- list(toPlot, lr.anova, p.vals.1, p.vals.2, cnt)
}
