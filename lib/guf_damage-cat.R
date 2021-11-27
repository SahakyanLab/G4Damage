

generateCat <- function(damage.type = c('PP', 'CPD', 'oxoG', 'cisplatin',
                                        'sonication', 'enzymatic', 'ancient'),
                        separator = ".", first = "nt"){
  
  rtn <- c()
  
  for (i in damage.type){
    if (i == 'cisplatin'){
      nt <- c('GG')
    } else if (i == 'PP'){
      nt <- c('TC', 'TT')
    } else if (i == 'CPD'){
      nt <- c('CT', 'TT')
    } else if (i == 'oxoG'){
      nt <- "G"
    } else if (i %in% c('sonication', 'ancient')) {
      nt <- "NN"
    } else if (i == 'enzymatic') {
      nt <- "NN"
    }
  
    cat <- rep("", length(nt))
    counter <- 1
    for (j in nt){
      if (first == "nt"){
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
