## FUNCTION ####################################################################
# This function takes a DNA sequence and looks for a regular expression.       #
#                                                                              #
# seq - a string or a vector of characters                                     #
# sequence argument (seq), the user can provide both a split DNA sequence (a   #
# vector of bases like c("A","T","G","G",...) ) and the merged one ("ATGG...").#
# The allowed characters are N, A, T, G, and C. Small letters are also allowed,#
# although be careful to account those in custom-defined x.strand.regex argu-  #
# ments.                                                                       #
#                                                                              #
# type - "GPQS"=="GQS", "APQS"=="AQS", "other"                                 #
# if type is "GPQS", the program will look for a classical G-quadruplex defini-#
# tion, with the maximum loop size additionally defined by the maxloopsize     #
# keyword. If "APQS", the adenine eqivalent will be considered. If "other", the#
# function will expect a non NULL plus.strand.regex and minus.strand.regex ar- #
# guments with the regular expressions to be searched in both + and - strands. #
#                                                                              #
# maxloopsize - numeric value                                                  #
# sets the maximum loop-size of the classical PQS motifs, used only when type  #
# is set to "GPQS" or "APQS".                                                  #
#                                                                              #
# plus.strand.regex - a string of regular expression                           #
# sets the regular expression to be searched in the sequence. The customly de- #
# fined values are read only if type is set to "other".                        #
#                                                                              #
# minus.strand.regex - a string of regular expression                          #
# sets the regular expression to be searched in the sequence. This one however #
# is the reverse and complementar equivalent of plus.strand.regex and is       #
# required for the program to account both strands. The customly defined values#
# are read only if type is set to "other".                                     #
#                                                                              #
# ... - aditional arguments are passed to the gregexp() function!              #
################################################################################
DNAregexfinder <- function(seq=seq, type="GPQS", maxloopsize=12,
                           plus.strand.regex=NULL, minus.strand.regex=NULL, ...){

  if(length(seq)>1){ # the string is split
    seq <- paste(seq, collapse="")
  }
  if(type=="GPQS" | type=="GQS"){
    plus.strand.regex  <- paste("([gG]{3}[NATGCnatgc]{1,",maxloopsize,"}){3}[gG]{3}",sep="")
    minus.strand.regex <- paste("([cC]{3}[NATGCnatgc]{1,",maxloopsize,"}){3}[cC]{3}",sep="")
  }
  if(type=="APQS" | type=="AQS"){
    plus.strand.regex  <- paste("([aA]{3}[NATGCnatgc]{1,",maxloopsize,"}){3}[aA]{3}",sep="")
    minus.strand.regex <- paste("([tT]{3}[NATGCnatgc]{1,",maxloopsize,"}){3}[tT]{3}",sep="")
  }
  if(type!="GPQS" & type!="APQS" & type!="GQS" & type!="AQS" & type!="other"){
    stop("QuadParser: non-recognisable type input.")
  }
  if(type=="other" & length(c(plus.strand.regex, minus.strand.regex))!=2){
    stop("QuadParser: type is set to 'other' but the regexes for both + and - strands are not specified.")  
  }

  plus.strand  <- gregexpr(text=seq, pattern=plus.strand.regex, ...)
  minus.strand <- gregexpr(text=seq, pattern=minus.strand.regex, ...)

  start.pos  <- c( as.vector(plus.strand[[1]]), as.vector(minus.strand[[1]]) )
  seq.length <- c( attr(plus.strand[[1]], "match.length"), attr(minus.strand[[1]], "match.length") )
  strand     <- c( rep("+",length(plus.strand[[1]])), rep("-",length(minus.strand[[1]])) )

  # Checking whether there are 0 returns (-1 by gregexpr)
  rm.ind <- which(start.pos==-1)
  if(length(rm.ind)!=0){
    start.pos  <- start.pos[-rm.ind]
    seq.length <- seq.length[-rm.ind]
    strand     <- strand[-rm.ind]
  }

  if(length(start.pos)!=0){ # there ARE detected occurrences
  
    # sorting the results in the order of their start.pos:
    new.order  <- order(start.pos)
  
    start.pos  <- start.pos[new.order]
    seq.length <- seq.length[new.order]
    strand     <- strand[new.order]   
    num.occ    <- length(start.pos)
    num.occ.plus  <- length(which(strand=="+"))
    num.occ.minus <- length(which(strand=="-"))
    sequence  <- sapply(1:num.occ, FUN=function(i){
                        substr(seq, start=start.pos[i], stop=(start.pos[i]+seq.length[i]-1) )
                       }, simplify=TRUE, USE.NAMES=FALSE)
  
  } else { # there are NO detected occurrences
    start.pos     <- 0
    seq.length    <- 0
    strand        <- 0   
    num.occ       <- 0
    num.occ.plus  <- 0
    num.occ.minus <- 0
    sequence      <- 0
  }

  QP.RESULTS <- NULL
  QP.RESULTS$start.pos     <- start.pos
  QP.RESULTS$seq.length    <- seq.length
  QP.RESULTS$strand        <- strand   
  QP.RESULTS$sequence      <- sequence
  QP.RESULTS$num.occ       <- num.occ
  QP.RESULTS$num.occ.plus  <- num.occ.plus
  QP.RESULTS$num.occ.minus <- num.occ.minus

  return(QP.RESULTS)

}
## FUNCTION ####################################################################
