################################################################################
# This function takes the starting/ending coordinates and a stratifying factor #
# (space) for the query segments and the subject segments. It then returns a   #
# matrix with all the indices of the overlapping segments.                     #
################################################################################
WhichOverlap <- function(start.query=start, end.query=end, space.query=chr,
                         start.subject=anno$txStart, end.subject=anno$txEnd,
                         space.subject=anno$chrom, maxgap=-1L, minoverlap=1L){

  #> source("http://bioconductor.org/biocLite.R")
  #> biocLite("IRanges")
  library(GenomicRanges)

  # RangedData function generates the corresponding objects, where it also reor-
  # ders the data, hence we provide subject.ind and query.ind index object in
  # order to account for the original indices for the data at the supplied order.
  
  # [ADIB] RangedData is deprecated in favor of GRanges() from GenomicRanges package.
  # subject <- RangedData( IRanges(start = start.subject, end = end.subject),
  #                        space = space.subject, subject.ind = 1:length(start.subject) )
  # 
  # query <- RangedData( IRanges(start = start.query, end = end.query),
  #                      space = space.query,  query.ind = 1:length(start.query) )

  space.subject <- unlist(strsplit(space.subject, "\\."))
  seqname.subject <- space.subject[seq(1, length(space.subject), 2)]
  strand.subject <- space.subject[seq(2, length(space.subject), 2)]
  
  subject <- GRanges(seqnames = seqname.subject,
                     ranges = IRanges(start = start.subject, end = end.subject),
                     strand = strand.subject)
  
  space.query <- unlist(strsplit(space.query, "\\."))
  seqname.query <- space.query[seq(1, length(space.query), 2)]
  strand.query <- space.query[seq(2, length(space.query), 2)]
  
  query <- GRanges(seqnames = seqname.query,
                     ranges = IRanges(start = start.query, end = end.query),
                     strand = strand.query)
  
  ol <- as.matrix(findOverlaps(query=query, subject=subject, type="any",
                            maxgap=maxgap, minoverlap=minoverlap, select="all"))

  # return( cbind(query=query[ol[,"queryHits"], ]$query.ind,
  #               subject=subject[ol[,"subjectHits"], ]$subject.ind) )

  colnames(ol) <- c("query", "subject")
  return(ol)
}
################################################################################
