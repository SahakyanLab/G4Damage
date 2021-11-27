
#cols = c(hue_pal()(3), hue_pal(l = 80, h.start = 180)(3), 'black')

colorChooser <- function(cats){
  
  rtn <- c()
  
  cols <- c(RColorBrewer::brewer.pal(7, "Set2"), "#864879", "#E9A6A6")
  
  for (cat in cats){
    
    if (grepl("CT.CPD", cat)) i <- 1
    else if (grepl("TT.CPD", cat)) i <- 2
    else if (grepl("TC.PP", cat)) i <- 3
    else if (grepl("TT.PP", cat)) i <- 4
    else if (grepl("GG.cisplatin", cat)) i <- 5
    else if (grepl("G.oxoG", cat)) i <- 6
    else if (grepl("NN.sonication", cat)) i <- 7
    else if (grepl("NN.ancient", cat)) i <- 8
    else if (grepl("NN.enzymatic", cat)) i <- 9
    
    rtn <- c(rtn, cols[i])
  }
  return(rtn)
}
