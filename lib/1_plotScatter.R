plotScatter <- function(altDNACode, altDNAStructure,
                        damage.types = c('PP', 'CPD'),
                        damage.strands = c('sense', 'antisense'),
                        dmg.at.altDNA.path,
                        ymax = 50,
                        writePDF = TRUE,
                        writePDF.height = 6,
                        writePDF.width = 8,
                        writePDF.path,
                        ...){


  # Dependencies: data.table, "GEN_scatterPlotInfo.R"

  dmg.at.altDNA = fread(dmg.at.altDNA.path)

  # [ADIB] What is this?
  dstrand.type <- function(s){
    if (s == "antisense"){
      return(3)
    }
    else if (s == "sense"){
      return(1)
    }
  }

  if (writePDF == TRUE){
    pdf(writePDF.path, height = writePDF.height, width = writePDF.width)
  }

  for (damage.type in damage.types){

    cats = generateCat(damage.type = damage.type)
    pdimers <- substr(cats, 1, 2)

    par(mar=c(5.1, 7.1, 4.1, 10.1), xpd=TRUE)

    plot(NULL,
         xlim = c(0,length(bins)), ylim = c(-ymax, ymax),
         main = paste(damage.type, " Damage in G-Quadruplex Sites", sep = ""),
         ylab = paste0("UV Damage at ", altDNAStructure,
                       " Sites \n (Percentage Change)"),
         xlab = paste0("Quality Score of ", altDNAStructure),
         xaxt = 'n') #, ...)


    axis(side = 1, at=0:length(bins),
         labels = paste0(as.character(c(0,bins)), "-",
                         as.character(c(bins, 100)), "%"))

    xs <- 0:length(bins)

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

      }

    }

    legend("topright", inset=c(-0.38,0),
           legend = c(paste0(substr(cats, 1,2), "-",
                             substr(cats, 4, nchar(cats))),
                      '' ,'Sense', 'Anti-sense'),
           lty = c(1,1,1,1,2),
           col = c(colorChooser(cats), 'white', 'black', 'black'))

  }

  while(!(names(dev.cur()) %in% c('null device', 'RStudioGD'))){
    dev.off()
  }

}
