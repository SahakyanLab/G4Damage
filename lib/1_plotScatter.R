plotScatter <- function(dmg.names = c('PP', 'CPD'),
                        dmg.strands = c('sense', 'antisense'),
                        dmg.at.G4.path,
                        bins,
                        raw.counts = FALSE,
                        cal.p.vals = FALSE,
                        ymax = 50,
                        combine.plot = FALSE,
                        writePDF.height = 6,
                        writePDF.width = 8,
                        writePDF.path = NULL){


  # Dependencies: data.table, "GEN_scatterPlotInfo.R"

  # if (!missing(dmg.at.G4.path)) {
  #   dmg.at.G4 <- fread(dmg.at.G4.path)
  # } else {
  #   dmg.at.G4 <- dmg.at.G4.dt
  # }
  
  dmg.at.G4 <- fread(dmg.at.G4.path)

  dstrand.type <- function(s){
    if (s == "antisense"){
      return(3)
    }
    else if (s == "sense"){
      return(1)
    }
  }

  if (!is.null(writePDF.path)){
    pdf(writePDF.path, height = writePDF.height, width = writePDF.width)
  }

  point.data <- "processed_data/points.tsv"
  #unlink(point.data)

  legend.lab <- NULL
  for (dmg.name in dmg.names){

    cats = generateCat(dmg.name = dmg.name)
    pdimers <- substr(cats, 1, 2)
    pdimers <- gsub("\\.$", "", pdimers)

    if (!combine.plot | dmg.name == dmg.names[1]) {

      par(mar=c(5.1, 7.1, 4.1, 10.1), xpd=TRUE)

      plot(NULL,
           xlim = c(0, length(bins)), ylim = c(-ymax, ymax),
           main = paste(dmg.name, " Damage in G-Quadruplex Sites", sep = ""),
           ylab = paste0("Damage propensity change (%)"),
           xlab = "G4 stability score (mm%)",
           xaxt = 'n', cex.axis = 0.8) #, ...)

      xs <- 0:length(bins)

      axis.lab <- paste0("(", as.character(c(0, bins)), ", ",
                         as.character(c(bins, 100)), "]")
      axis.lab[1] <- sub("\\(", "[", axis.lab[1])
      axis(side = 1, at=0:length(bins), labels = NA)
      text(x = xs,
           y = par("usr")[3] - 0.08 * (par("usr")[4] - par("usr")[3]),
           labels = axis.lab,
           xpd = NA,
           srt = 35,
           cex = 0.8)
    }

    sig.plots <- NULL
    for (pdimer in pdimers){

      for (dmg.strand in dmg.strands){
        
        cat(dmg.name, "-", pdimer, " at ", dmg.strand, " stand\n", sep = "")
        
        toPlot <- scatterPlotInfo(dmg.at.G4 = dmg.at.G4,
                                  strand = dmg.strand,
                                  dmg.name = dmg.name,
                                  pdimer = pdimer,
                                  bins = bins,
                                  raw.counts = raw.counts,
                                  cal.p.vals = cal.p.vals)
        cat("\n")
        sig.plots[[paste0(pdimer, "-", dmg.strand)]] <- toPlot
        
        print(toPlot[[1]])
        
        points(xs, toPlot[[1]],
               col = colorChooser(paste(pdimer, dmg.name, sep=".")),
               type = "l", lty = dstrand.type(dmg.strand), lwd = 3)
        points(xs, toPlot[[1]],
               col = colorChooser(paste(pdimer, dmg.name, sep=".")),
               pch = 19)

        fwrite(data.table(damage = paste0(dmg.name, "-", pdimer),
                          strand = dmg.strand,
                          bin = xs,
                          Percent = toPlot[[1]]),
               point.data, sep = "\t", append = TRUE)

      }

    }

    if (cal.p.vals) {
      sig.codes <- c("***" = 0.001, "**" = 0.01, "*" = 0.05, " " = 1)
      for (pdimer in pdimers){
        for (dmg.strand in dmg.strands){
          #a <<- sig.plots[[paste0(pdimer, "-", dmg.strand)]]
          text(x = xs, y = sig.plots[[paste0(pdimer, "-", dmg.strand)]][[1]],
               labels = sapply(sig.plots[[paste0(pdimer, "-", dmg.strand)]][[3]],
                               function(p) {
                                 s <- sig.codes[p <= sig.codes] |> names()
                                 s[1]
                               }))
          text(x = xs, y = sig.plots[[paste0(pdimer, "-", dmg.strand)]][[1]],
               labels = sapply(sig.plots[[paste0(pdimer, "-", dmg.strand)]][[4]],
                               function(p) {
                                 s <- sig.codes[p <= sig.codes] |> names()
                                 s[1]
                               }),
               #pos = 2,
               adj = c(1.5, NA),
               col = "grey50")
          text(x = xs[length(xs)],
               y = sig.plots[[paste0(pdimer, "-", dmg.strand)]][[1]][length(xs)],
               labels = names(sig.codes[
                 sig.plots[[paste0(pdimer, "-", dmg.strand)]][[2]] <= sig.codes])[1],
               #pos = 4,
               adj = c(-0.5, NA),
               col = "grey")
        }
      }
    }
    

    if (combine.plot) {
      lab <- paste0(pdimers, "-", strsplit(cats, "\\.")[[1]][2])
      legend.lab[(((l <- length(legend.lab))+1):(l+length(lab)))] <- lab
    } else {
      legend.lab <- paste0(pdimers, "-", strsplit(cats, "\\.")[[1]][2])
    }

    if (!combine.plot | dmg.name == dmg.names[length(dmg.names)]) {
      legend("topright", inset=c(-0.38,0),
             legend = c(legend.lab,
                        if(length(dmg.strands) == 2) c('' ,'Sense', 'Anti-sense')),
             lty = c(rep(1, length(legend.lab)), if(length(dmg.strands) == 2) c(1, 1,2)),
             lwd = 1.3,
             seg.len = 4,
             col = c(colorChooser(gsub("-", ".", legend.lab)),
                     if(length(dmg.strands) == 2) c('white', 'black', 'black')),
             cex = 0.8)
    }

  }

  while(!(names(dev.cur()) %in% c('null device', 'RStudioGD'))){
    dev.off()
  }

}
