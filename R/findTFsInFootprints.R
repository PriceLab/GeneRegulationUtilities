# findTFsInFootprints.R
#------------------------------------------------------------------------------------------------------------------------
library(GenomicRanges)
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
findTFsInFootprints <- function(tbl.fp, tbl.motifToGenes, chromosome, startLoc, endLoc, transcriptionFactors=NA)
{
    if(all(is.na(transcriptionFactors)))
        transcriptionFactors <- sort(unique(tbl.motifToGenes$tfs))

    stopifnot(all(c("chr", "mfpStart", "mfpEnd", "motifName") %in% colnames(tbl.fp)))

    tbl.sub <- subset(tbl.fp, chr==chromosome & mfpStart >= startLoc & mfpEnd <= endLoc)
    motifs.per.fp <- tbl.sub$motifName
    tbl.hits <- subset(tbl.motifToGenes, motif %in% motifs.per.fp)
    indicesMatched <- match(transcriptionFactors, tbl.hits$tfs)
    tbl.hitAndMatched <- subset(tbl.hits, tfs %in% transcriptionFactors)
    motifs <- unique(tbl.hitAndMatched$motif)
    tbl.sub <- subset(tbl.sub, motifName %in% motifs)

    tfs.per.motif <- lapply(motifs, function(m) subset(tbl.hitAndMatched, motif==m)$tfs)
    names(tfs.per.motif) <- motifs
    tfsPerRow <- unlist(lapply(tbl.sub$motifName, function(motifName)
        paste(intersect(transcriptionFactors, subset(tbl.motifToGenes, motif==motifName)$tfs), collapse=";")))
    tbl.augmented <- cbind(tbl.sub, tfsMatched=tfsPerRow, stringsAsFactors=FALSE)

    invisible(tbl.augmented)

} # findTFsInFootprints
#------------------------------------------------------------------------------------------------------------------------
# identify some proximal upstream snps, within
# with mef2c at hg39 chr5:88,716,241-88,906,105, coded on the minus strand, here are some proximal upstream snps
# reported by mariette
#
# true overlaps with mfpStart-mfpEnd: rs3814428, rs770463, rs446500
# near misses: rs10038371 rs214136 rs16876775
#
# fivenum(tbl$GRCh38.liftover) [1] 88640561 88740812 88828020 88921378 89001975
#   dim(subset(tbl, GRCh38.liftover > 88906105))  [1] 43 18
#   tbl.upstream <- tbl[which(tbl$GRCh38.liftover - 88906105 > 0),]   # 43 18
#   upstream.loc <- tbl.upstream$GRCh38.liftover - 88906105
# head(tbl.upstream[, c(1,3,4,5,6,17:19)], n=20)
#     Symbol CHR        SNP       BP A1 GRCh38.liftover   GRCh37 upstream
# 115  MEF2C   5 rs10085009 88241752  A        88910179 88205996     4074
# 150  MEF2C   5 rs13357047 88246411  A        88914838 88210655     8733
# 119  MEF2C   5   rs770466 88251210  C        88919637 88215454    13532
# 81   MEF2C   5   rs304132 88251350  A        88919777 88215594    13672
# 153  MEF2C   5 rs10044186 88254553  C        88922980 88218797    16875
# 26   MEF2C   5 rs10037047 88255720  G        88924147 88219964    18042
# 44   MEF2C   5  rs190982* 88259176  G        88927603 88223420    21498
# 42   MEF2C   5   rs301717 88267164  A        88935591 88231408    29486
# 30   MEF2C   5   rs301718 88267595  G        88936022 88231839    29917
# 34   MEF2C   5   rs301719 88268671  C        88937098 88232915    30993
# 40   MEF2C   5   rs301720 88269410  A        88937837 88233654    31732
# 35   MEF2C   5   rs301723 88271430  A        88939857 88235674    33752
# 36   MEF2C   5   rs301725 88272216  T        88940643 88236460    34538
# 29   MEF2C   5  rs7700950 88272698  T        88941125 88236942    35020
# 33   MEF2C   5   rs301726 88272826  T        88941253 88237070    35148
# 31   MEF2C   5   rs301727 88273286  G        88941713 88237530    35608
# 32   MEF2C   5   rs301728 88273557  T        88941984 88237801    35879
# 39   MEF2C   5   rs301729 88273797  G        88942224 88238041    36119
# 41   MEF2C   5  rs1650626 88274962  G        88943389 88239206    37284
# 37   MEF2C   5  rs1650627 88275084  G        88943511 88239328    37406
findSNPsInFootprints <- function(tbl.fp, tbl.motifToGenes,
                                 chromosome, region.startLoc, region.endLoc,
                                 snp.chromosome, snp.loc,
                                 transcriptionFactors=NA,
                                 padding=0)
{
    tbl.fpSub <- findTFsInFootprints(tbl.fp, tbl.motifToGenes,
                                     chromosome, region.startLoc, region.endLoc,
                                     transcriptionFactors)
    if(nrow(tbl.fpSub) == 0)
        return(data.frame())

    gr.fp <- with(tbl.fpSub, GRanges(seqnames=chr, IRanges(start=mfpStart, end=mfpEnd)))
    gr.snp <- GRanges(seqnames=snp.chromosome, IRanges(start=snp.loc-padding, end=snp.loc+padding))
    #browser()
    tbl.overlaps <- as.data.frame(findOverlaps(gr.fp, gr.snp))
    snp <- rep(-1, length(gr.fp))
    snp[tbl.overlaps$queryHits] <- snp.loc[tbl.overlaps$subjectHits]

    tbl.fpSub$snp <- snp
    result <- subset(tbl.fpSub, snp != -1)
    #browser()
    result[, c("chr", "mfpStart", "mfpEnd", "motifStart", "motifEnd", "sequence", "motifName", "snp", "tfsMatched")]

} # findTFsInFootprints
#------------------------------------------------------------------------------------------------------------------------
# me2fc upstream mariette-reported snps in footprints
#    seqnames    start      end width strand motifName score seqnames.1  start.1    end.1 width.1 strand.1
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

