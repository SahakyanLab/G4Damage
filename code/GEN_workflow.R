## ---- Description


## ---- setWD

#setwd("~/Actual Work/WIMM/g-quadruplex paper/")

#on ox network:
#setwd("~/CLAUDIA_WORKSPACE/g4_paper/")

#on ox computer:
setwd("~/Documents/g4_paper/")

## ---- setParameters
#parameters that will likely need to be changed

NCPU <- 1

altDNACode <- "g4"#ox8
altDNAStructure <- "G-Quadruplex"#"8-Oxoguanine"

highLow.cutoff <- 19
bins <- c(5, 10, 15, 20, 25, 30, 35, 40)


## ---- do

toDo <- new.env()

toDo$generateTable <- TRUE

toDo$countUVDamages <- TRUE

toDo$scatterPlot <- TRUE
savePDF.scatter = TRUE

toDo$categorizeUVDamages <- TRUE

## ---- fixedParams

chromosome.list <- c(1:22, 'X', 'Y')#ignore mitochondrial dna

damage.strands <- c('antisense', 'sense')
damage.types <- c('CPD', "PP")

## ---- dependentFunctions

source("code/GEN_findOligonucleotideCounts.R")
source("code/GEN_countDamageOverlap.R")

source("code/GEN_damage-cat.R")
source("code/GEN_experiment-info.R")
source("code/GEN_smooth.R")
source("code/GEN_colorChooser.R")

#dependent functions from the TrantorR library

TrantoRLib <- "code/TrantoRextr/"
source(paste(TrantoRLib, "GEN_readfasta.R", sep=""))
source(paste(TrantoRLib, "UTIL_readLinesFast.R", sep=""))
source(paste(TrantoRLib, "GEN_loadGenome.R", sep=""))
source(paste(TrantoRLib, "GEN_get_revcomp.R", sep=""))
source(paste(TrantoRLib, "GEN_WhichOverlap.R", sep=""))
rm(TrantoRLib)

#also need libraries
library(Biostrings)
library(scales)


## ---- parallelize

library(doParallel)
source(file = "code/GEN_doParallelsetup.R")



## ---- startCode

#####################################################################################################################################
############################################################   start code   #########################################################
#####################################################################################################################################

## ---- generateTable
## generate table with raw G4/OX8 data that counts the number of potential UV damage sites there are

##data table should have the following columns:
#chr: chromosome where the alternative DNA structure was observed
#genomic.start: genomic start coordinate
#genomic.end: genomic end coordinate
#relseq: sequence on the strand where the alternative DNA structure was observed
#mm: quality score
#also any additional columns are fine (can be in any order)

if (toDo$generateTable == TRUE){
  source(file = paste("code/", toupper(altDNACode), "_generate", toupper(altDNACode), "Table.R", sep=""))
  eval(parse(text = paste("generate", toupper(altDNACode), "Table()", sep="")))
  rm(list = ls(pattern = paste("generate", toupper(altDNACode), "Table", sep="")))
}



## ---- UVDamageatAltDNASites

if (toDo$countUVDamages == TRUE){
  source(file = "code/GEN_UVDamageAtAltDNACounts.R")
  UVDamageAtaltDNACounts(altDNACode = altDNACode)
  rm(UVDamageAtaltDNACounts)
}


## ---- scatterPlots

if (toDo$scatterPlot == TRUE){
  source("code/GEN_plotScatter.R")
  plotScatter(altDNACode = altDNACode,
              writePDF = savePDF.scatter)
  rm(savePDF.scatter)
}






