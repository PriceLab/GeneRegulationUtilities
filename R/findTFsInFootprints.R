# findTFsInFootprints.R
#------------------------------------------------------------------------------------------------------------------------
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
#if(!exists("tbl.fpAnnotated")){
#   printf("loading serialized version of uniq_fimo_fp.bed");
#   elapsed.time <- system.time(load(system.file(package="GeneRegulationUtilities", "extdata", "tbl.fpAnnotated.RData")))[["elapsed"]]
#   printf("time to load %s: %s seconds", "tbl.fpAnnotated", elapsed.time)
#   }
#if(!exists("tbl.motifToMultipleGenes")){
#   printf("loading serialized version of motif_to_tf_mappings_with_tfclass_include_multiple.csv");
#   elapsed.time <- system.time(load(system.file(package="GeneRegulationUtilities", "extdata","tbl.motifToMultipleGenes.RData")))[["elapsed"]]
#   printf("time to load %s: %s seconds", "tbl.motifToMultipleGenes", elapsed.time)
#   }
#
#------------------------------------------------------------------------------------------------------------------------
findTFsInFootprints <- function(tbl.fp, tbl.motifToGenes, chromosome, startLoc, endLoc, transcriptionFactors)
{
    tbl.sub <- subset(tbl.fp, chr==chromosome & start >= startLoc & end <= endLoc)
    motifs.per.fp <- tbl.sub$name
    tbl.hits <- subset(tbl.motifToGenes, motif %in% motifs.per.fp)
    tbl.hitAndMatched <- subset(tbl.hits, tfs %in% transcriptionFactors)
    motifs <- unique(tbl.hitAndMatched$motif)
    tbl.sub <- subset(tbl.sub, name %in% motifs)

    tfs.per.motif <- lapply(motifs, function(m) subset(tbl.hitAndMatched, motif==m)$tfs)
    names(tfs.per.motif) <- motifs
    tfsPerRow <- unlist(lapply(tbl.sub$name, function(motifName)
        paste(intersect(transcriptionFactors, subset(tbl.motifToGenes, motif==motifName)$tfs), collapse=";")))
    tbl.augmented <- cbind(tbl.sub, tfsMatched=tfsPerRow, stringsAsFactors=FALSE)

    invisible(tbl.augmented)

} # findTFsInFootprints
#------------------------------------------------------------------------------------------------------------------------
# me2fc upstream mariette-reported snps in footprints
#    seqnames    start      end width strand      name score seqnames.1  start.1    end.1 width.1 strand.1
# 1      chr5 88669200 88669211    12      -  MA0486.2  50.7       chr5 88669204 88669204       1        *
# 2      chr5 88669200 88669211    12      -  MA0771.1  50.1       chr5 88669204 88669204       1        *
# 3      chr5 88884580 88884599    20      -  MA0528.1  56.3       chr5 88884588 88884588       1        *
# 4      chr5 88884581 88884590    10      + KLF16_DBD  51.1       chr5 88884588 88884588       1        *
# 5      chr5 88884581 88884590    10      +  MA0079.3  60.5       chr5 88884588 88884588       1        *
# 6      chr5 88884581 88884594    14      +  MA0516.1  51.6       chr5 88884588 88884588       1        *
# 7      chr5 88884581 88884590    10      +  MA0741.1  51.1       chr5 88884588 88884588       1        *
# 8      chr5 88884581 88884590    10      +  MA0746.1  50.0       chr5 88884588 88884588       1        *
# 9      chr5 88884583 88884589     7      -    MAZ.p2  50.9       chr5 88884588 88884588       1        *
# 10     chr5 88884584 88884591     8      -  GTF2I.p2  50.2       chr5 88884588 88884588       1        *
# 11     chr5 88884585 88884593     9      +  MA0753.1  50.1       chr5 88884588 88884588       1        *
# 12     chr5 88900174 88900187    14      -  MA0505.1  53.7       chr5 88900187 88900187       1        *
#------------------------------------------------------------------------------------------------------------------------
prepareData <- function()
{
   f <- 'tfdb_1mb_GRCh38.83.bed'  # too big
   f <- 'tiny.bed'                # for testing
   f <- 'uniq_fimo_fp.bed'        # using this one

   tbl.fpAnnotated <- read.table(f, header=FALSE, sep="\t", as.is=TRUE)
   new.colnames <- c("chr", "start", "end", "score", "strand", "name", "info", "sequence", "width", "chrom", "bigStart",
                     "bigEnd", "x", "y", "z", "zz")
   colnames(tbl.fpAnnotated) <- new.colnames
   tbl.fpAnnotated$chr <- paste("chr", tbl.fpAnnotated$chr, sep="")
   save(tbl.fpAnnotated, file = "tbl.fpAnnotated.RData")

  tbl.motifToMultipleGenes <- read.table("motif_to_tf_mappings_with_tfclass_include_multiple.csv", sep=",", header=TRUE, as.is=TRUE)
  save(tbl.motifToMultipleGenes, file="tbl.motifToMultipleGenes.RData")


} # prepareData
#------------------------------------------------------------------------------------------------------------------------

