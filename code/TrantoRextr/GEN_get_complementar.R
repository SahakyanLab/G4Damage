###############################################################################
get.complementar <- function(seq="GCATTCCACA"){
  seq <- strsplit(seq,"")[[1]]
  seq2 <- seq
  seq2[which(seq=="A")] <- "T"
  seq2[which(seq=="G")] <- "C"
  seq2[which(seq=="T")] <- "A"
  seq2[which(seq=="C")] <- "G"
  return(paste(seq2,collapse=""))
}
###############################################################################

