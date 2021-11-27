plotScatter <- function(altDNACode, altDNAStructure,
                        damage.types = c('PP', 'CPD'),
                        damage.strands = c('sense', 'antisense'),
                        dmg.at.altDNA.path,
                        ymax = 50,
                        combine.plot = FALSE,
                        writePDF.height = 6,
                        writePDF.width = 8,
                        writePDF.path = NULL,
                        ...){


  # Dependencies: data.table, "GEN_scatterPlotInfo.R"

  dmg.at.altDNA = fread(dmg.at.altDNA.path)

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
  for (damage.type in damage.types){

    cats = generateCat(damage.type = damage.type)
    pdimers <- substr(cats, 1, 2)
    pdimers <- gsub("\\.$", "", pdimers)

    if (!combine.plot | damage.type == damage.types[1]) {
      par(mar=c(5.1, 7.1, 4.1, 10.1), xpd=TRUE)
      
      plot(NULL,
           xlim = c(0,length(bins)), ylim = c(-ymax, ymax),
           main = paste(damage.type, " Damage in G-Quadruplex Sites", sep = ""),
           ylab = paste0("Damage propensity change (%)"),
           xlab = paste0(altDNACode, " stability score (mm%)"),
           xaxt = 'n', cex.axis = 0.8) #, ...)
      
      xs <- 0:length(bins)
      
      axis.lab <- paste0("(", as.character(c(0,bins)), ", ",
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

    for (pdimer in pdimers){

      for (damage.strand in damage.strands){

        toPlot <- scatterPlotInfo(altDNACode = altDNACode,
                                  dmg.at.altDNA = dmg.at.altDNA,
                                  strand = damage.strand,
                                  damage.type = damage.type, pdimer = pdimer,
                                  bins = bins) #, ...)
        
        points(xs, toPlot,
               col = colorChooser(paste(pdimer, damage.type, sep=".")),
               type="l", lty = dstrand.type(damage.strand), lwd = 3)
        points(xs, toPlot,
               col = colorChooser(paste(pdimer, damage.type, sep=".")),
               pch = 19)
        
        fwrite(data.table(damage = paste0(damage.type, "-", pdimer),
                          strand = damage.strand,
                          bin = xs,
                          Percent = toPlot),
               point.data, sep = "\t", append = TRUE)

      }

    }
    
    if (combine.plot) {
      legend.lab[damage.type] <- paste0(pdimers, "-", strsplit(cats, "\\.")[[1]][2])
    } else {
      legend.lab <- paste0(pdimers, "-", strsplit(cats, "\\.")[[1]][2])
    }

    if (!combine.plot | damage.type == damage.types[length(damage.types)]) {
      legend("topright", inset=c(-0.38,0),
             legend = c(legend.lab,
                        if(length(damage.strands) == 2) c('' ,'Sense', 'Anti-sense')),
             lty = c(rep(1, length(legend.lab)), if(length(damage.strands) == 2) c(1, 1,2)),
             lwd = 2,
             col = c(colorChooser(gsub("-", ".", legend.lab)),
                     if(length(damage.strands) == 2) c('white', 'black', 'black')),
             cex = 0.8)
    }

  }

  while(!(names(dev.cur()) %in% c('null device', 'RStudioGD'))){
    dev.off()
  }

}
