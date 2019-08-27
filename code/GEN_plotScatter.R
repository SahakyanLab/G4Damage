plotScatter <- function(altDNACode,
                        damage.types = c('PP', 'CPD'), damage.strands = c('sense', 'antisense'),
                        uv.at.altDNA.path = paste0("data/", altDNACode, "/", "uv-damage-at-", altDNACode, "-counts.csv"),
                        ymax = 50,
                        writePDF = TRUE, writePDF.height = 6, writePDF.width = 8,
                        writePDF.path = paste0("figures/", toupper(altDNACode), "_uv-damage-at-", altDNACode, "-sites.pdf"),
                        ...){
  
  
  source("code/GEN_scatterPlotInfo.R")
  
  uv.at.altDNA = data.table::fread(uv.at.altDNA.path)
  
  dstrand.type <- function(s){if (s == "antisense"){return(3)} else if (s == "sense"){return(1)}}
  
  if (writePDF == TRUE){pdf(writePDF.path, height = writePDF.height, width = writePDF.width)}
  
  for (damage.type in damage.types){
    
    cats = generateCat(ALL = FALSE, damage.type = damage.type)
    pdimers <- substr(cats, 1, 2)
   
    par(mar=c(5.1, 7.1, 4.1, 10.1), xpd=TRUE)
    
    plot(NULL,
         xlim = c(0,length(bins)), ylim = c(-ymax, ymax),
         main=paste(damage.type, " Damage in G-Quadruplex Sites", sep = ""), 
         ylab=paste0("UV Damage at ", altDNAStructure, " Sites \n (Percentage Change)"),
         xlab=paste0("Quality Score of ", altDNAStructure),
         xaxt = 'n', ...)
    
    
    axis(side = 1, at=0:length(bins), 
         labels = paste(as.character(c(0,bins)), "-", as.character(c(bins, 100)), "%", sep=""))
    
    xs <- 0:length(bins)
    
    for (pdimer in pdimers){
      
      for (damage.strand in damage.strands){
        
        toPlot <- scatterPlotInfo(altDNACode = altDNACode,
                                  uv.at.altDNA = uv.at.altDNA,
                                  strand = damage.strand,
                                  damage.type = damage.type, pdimer = pdimer,
                                  bins = bins, ...)
        points(xs, toPlot,
               col = colorChooser(paste(pdimer, damage.type, sep=".")), 
               type="l", lty = dstrand.type(damage.strand), lwd = 3)
        points(xs, toPlot,
               col = colorChooser(paste(pdimer, damage.type, sep=".")), pch = 19)
        
      }
      
    }
    
    legend("topright", inset=c(-0.35,0),
           legend = c(paste(substr(cats, 1,2), "-", substr(cats, 4, nchar(cats)), sep=""), '' ,'Sense', 'Anti-sense'),
           lty = c(1,1,1,1,2), col = c(colorChooser(cats), 'white', 'black', 'black'))
    
  }
  
  while(!(names(dev.cur()) %in% c('null device', 'RStudioGD'))){dev.off()}
  
}



