
#cols = c(hue_pal()(3), hue_pal(l = 80, h.start = 180)(3), 'black')

colorChooser <- function(cats){
  
  rtn <- c()

  cols <- c("#E77D72", "#53B74C", "#F2A49A", "#6BE079", "#6F9BF8", "#F4BF62", "#9C5A33", "#984EA3", "#171717")
  
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
