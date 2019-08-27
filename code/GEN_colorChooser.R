
#cols = c(hue_pal()(3), hue_pal(l = 80, h.start = 180)(3), 'black')

colorChooser <- function(cats, ls = c(65, 80),
                         cols = c(hue_pal(l = ls[1])(3), hue_pal(l = ls[2], h.start = 180)(3), 'black')){
  
  rtn <- c()
  
  for (cat in cats){
    
    if (cat == 'ALL'){
      i <- 7
    } else {
      if (nchar(cat) < 4 || nchar(cat) > 7){#if cat has less than 4 characters, it is either "CPD" or "PP"
        if (length(grep("CPD", cat))==0){#cat = PP
          i <- 3
        } else{
          i <- 6#cat = CPD
        }
      } else {
        if (length(grep("CPD", cat))==0){#cat = PP
          if (length(grep("TT", cat))==0){#cat = TT-PP
            i <- 1
          } else {
            i <- 2
          }
        } else{
          if (length(grep("TT", cat))==0){#cat = TT-PP
            i <- 4
          } else {
            i <- 5
          }
        }
      }
    }
    rtn <- c(rtn, cols[i])
  }
  return(rtn)
}
