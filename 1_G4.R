################################################################################
#
# Workflow to get similar plots in the publication
# 1. UV damage
#    dmg.names = c('CPD', 'PP')
#    dmg.type = "UV"
#    dmg.pattern = c("CT", "TC", "TT")
#    strand.sensitive = T
#    combine.plot = F
# 2. cisplatin damage
#    dmg.names = "cisplatin"
#    dmg.type = "cisplatin"
#    dmg.pattern = "GG"
#    strand.sensitive = T
# 3. 8oxoG damage
#    dmg.names = "oxoG"
#    dmg.type = "oxoG"
#    dmg.pattern = "G"
#    strand.sensitive = T
# 4. breakage damage
#    dmg.names = c("sonication", "enzymatic", "ancient")
#    dmg.type = "breakage"
#    dmg.pattern = "NN"
#    strand.sensitive = F
#    combine.plot = T
#
## Configuration ###############################################################
NCPU = 1

dmg.names = c("sonication", "enzymatic") #c("CPD", "PP") "oxoG" "cisplatin" c("sonication", "enzymatic") "ancient"
dmg.type = "breakage" # Selection: "cisplatin" "UV" "oxoG" "breakage"
dmg.pattern = "NN" # c("CT", "TC", "TT") "G" "GG" "NN"

strand.sensitive = T
dmg.strands = if (strand.sensitive) c("antisense", "sense") else "sense"

highLow.cutoff = 19
bins = c(5, 10, 15, 20, 25, 30, 35, 40)


## Task ########################################################################
toDo = NULL

toDo$generateTable         = F

toDo$countDamages          = F
include.weight = F # must have weight column

toDo$scatterPlot           = T
combine.plot = F
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
source("lib/1_DamageAtG4Counts.R")
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

  damageAtG4Counts(dmg.names = dmg.names,
                   G4.path = processed.PQSdata.path,
                   strand.sensitive = strand.sensitive,
                   include.weight = include.weight,
                   saveTab.path = dmg.cnt.G4.path)

}

#-------------------------------------------------------------------------------
if (toDo$scatterPlot){

  plotScatter(dmg.names = dmg.names,
              dmg.strands = dmg.strands,
              dmg.at.G4.path = dmg.cnt.G4.path,
              ymax = 80,
              combine.plot,
              writePDF.path = plot.path)
}
