################################################################################
## Configuration ###############################################################
NCPU = 1

altDNACode = "G4"
altDNAStructure = "G-Quadruplex"

dmg.type = "breakage" # Selection: cisplatin UV oxoG breakage
dmg.pattern = "NN" # c("CC", "CT", "TC", "TT") G GG NN 
strand.sensitive = F
include.weight = T # must have weight column

highLow.cutoff = 19
bins = c(5, 10, 15, 20, 25, 30, 35, 40)

damage.strands = if (strand.sensitive) c('antisense', 'sense') else 'sense'
damage.types = c('sonication', 'enzymatic', 'ancient') #c('CPD', 'PP') oxoG cisplatin c('sonication', 'enzymatic', 'ancient')
combine.plot = T

## Task ########################################################################
toDo = NULL

toDo$generateTable         = F

toDo$countDamages          = F

toDo$scatterPlot           = T
savePDF.scatter = T

## File path ###################################################################
PQSdata.path = "raw_data/PQSdata.txt"

# Outputs
processed.PQSdata.path = paste0("processed_data/", dmg.type, "_G4_data.csv")
dmg.cnt.G4.path = paste0("processed_data/", dmg.type, "-damage-at-G4-counts.csv")
plot.path = paste0("processed_data/", dmg.type, "-damage-at-G4-counts.pdf")


## Dependant functions #########################################################
source("lib/1_generateG4Table.R")
source("lib/1_findOligonucleotideCounts.R")
source("lib/1_DamageAtAltDNACounts.R")
source("lib/1_countDamageOverlap.R")
source("lib/1_plotScatter.R")
source("lib/1_scatterPlotInfo.R")

source("lib/guf_colorChooser.R")
source("lib/guf_damage-cat.R")

# Dependant functions from the TrantorR library
TrantoRLib <- "lib/TrantoRextr/"
source(paste(TrantoRLib, "GEN_WhichOverlap.R", sep=""))

## Dependant libraries #########################################################
suppressPackageStartupMessages( library(doParallel)    )
suppressPackageStartupMessages( library(Biostrings)    )
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(data.table)    )
suppressPackageStartupMessages( library(scales    )    )


#suppressPackageStartupMessages( library(venneuler)     )
#suppressPackageStartupMessages( library(stringi)       )

## Parallel setup ##############################################################
#source("lib/GEN_doParallelsetup.R")

## Directory setup #############################################################
suppressWarnings(dir.create("processed_data"))

## Main Code ###################################################################

#-------------------------------------------------------------------------------
if (toDo$generateTable){

  # To generate table with raw G4 data that counts the number of potential UV
  # or oxidative damage sites there are

  # The input data table should have the following columns:
  # chr           : chromosome where the alternative DNA structure was observed
  # genomic.start : genomic start coordinate
  # genomic.end   : genomic end coordinate
  # relseq        : sequence on the strand where the alternative DNA structure
  #                 was observed
  # mm            : quality score
  # also any additional columns are fine (can be in any order)
  #
  # The output data table will have the following columns:
  #   Chr, Quality, Strand, Start, End, RelSeq, RelSeq_pattern_cnt, RevComp,
  #   RevComp_pattern_cnt

  generateG4Table(G4.data.path = PQSdata.path,
                  dmg.pattern = dmg.pattern,
                  saveTab.path = processed.PQSdata.path)
}

#-------------------------------------------------------------------------------
if (toDo$countDamages){
  # count reported DNA damage on G4 structure.
  # Output table has the following columns:
  #  Chr, Strand, Start, End, Quality,      [G4 info]
  #  <dmg_type>.<dmg_nt>.sense.altDNA,      [count]
  #  <dmg_type>.<dmg_nt>.sense.overlap,     [count]
  #  <dmg_type>.<dmg_nt>.antisense.altDNA,  [count]
  #  <dmg_type>.<dmg_nt>.antisense.overlap  [count]
  #  ... (add more columns for more damage type and/or damage nucleotide)

  damageAtaltDNACounts(damage.types = damage.types,
                       altDNAData.path = processed.PQSdata.path,
                       strand.sensitive = strand.sensitive,
                       include.weight = include.weight,
                       saveTab.path = dmg.cnt.G4.path)

}

#-------------------------------------------------------------------------------
if (toDo$scatterPlot){

  plotScatter(altDNACode = altDNACode,
              altDNAStructure = altDNAStructure,
              damage.types = damage.types,
              damage.strands = damage.strands,
              dmg.at.altDNA.path = dmg.cnt.G4.path,
              ymax = 80,
              combine.plot,
              writePDF.path = plot.path)
}
