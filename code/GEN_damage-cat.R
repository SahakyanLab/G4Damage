

generateCat <- function(damage.type = c('PP', 'CPD'),
                        ALL=TRUE, separator = ".",
                        first = "pdimer"){
  
  rtn <- c()
  
  for (i in damage.type){
    if (ALL == TRUE){
      pdimers <- c('CC', 'CT', 'TC', 'TT')
    } else {
      if (i=='PP'){
        pdimers <- c('TC', 'TT')
      } else if (i == 'CPD'){
        pdimers <- c('CT', 'TT')
      } 
    }
    
    cat <- rep("", length(pdimers))
    counter <- 1
    for (j in pdimers){
      if (first == "pdimer"){
        cat[counter] <- paste(j, separator, i, sep="")
      } else if (first == "damage.type"){
        cat[counter] <- paste(i, separator, j, sep="")
      }
      counter <- counter+1
    }
    rtn <- c(rtn, cat)
  }
  
  return(rtn)
}
