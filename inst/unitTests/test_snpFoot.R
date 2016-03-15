library(RUnit)
library(snpFoot);
library(PrivateCoryData)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.fpAnnotated")){
   printf("loading tbl.fpAnnotated from PrivateCoryData package...")
   data(tbl.fpAnnotated)
   }

if(!exists("tbl.fpAnnotated")){
   printf("loading tbl.fpAnnotated from PrivateCoryData package...")
   data(tbl.fpAnnotated)
   }

if(!exists("tbl.motifToMultipleGenes")){
   printf("loading tbl.motifToMultipleGenes from PrivateCoryData package...")
   data(tbl.motifToMultipleGenes)
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_1_footprint_SAT2B_tf()
   test_1_footprint_two_tfs()
   test_emptyResults()
   test_360k_chr5_alzheimers_region()
   test_findSNPsInFootprints()
   test_findSNPsNearFootprints()
   test_LXH1_matching()
   #test_snpFootDisplay()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
# tbl.motifToMultipleGenes (from seth's "motif_to_tf_mappings_with_tfclass_include_multiple.csv" tells us that
# that SATB2 binds to pwm motif MA0679.1.
# strategy:
#  1) find a footprint in tbl.fpAnnotated by direct inspection: grep("MA0679", tbl.fpAnnotated$info)[1]  #  891
#  2) use the chromosomal location of that cherry-picked footprint
#  3) make sure that findTFsInFootprints retrieves that footprint, and annotates it to SATB2
test_1_footprint_SAT2B_tf <- function()
{
  printf("--- test_1_footprint_SAT2B_tf")

  loc <- tbl.fpAnnotated[grep("MA0679", tbl.fpAnnotated$motifName)[1], c("chr", "mfpStart", "mfpEnd")]
  target.tf <- "SATB2"
  tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                             chromosome=loc$chr, startLoc=loc$mfpStart, endLoc=loc$mfpEnd,
                             transcriptionFactors=target.tf)

    # tbl.fpAnnotated identifies 3 motifs in this one footprint, and each of them is associated with SATB2
    # we therefore expect 3 rows, each describing the same footprint, and each with a different motif

  checkEquals(dim(tbl), c(3, 18))
  checkEquals(tbl$tfsMatched, rep("SATB2", 3))
  checkEquals(tbl$motifName, c("MA0679.1", "MA0756.1", "MA0757.1"))

} # test_1_footprint_SAT2B_tf
#------------------------------------------------------------------------------------------------------------------------
test_1_footprint_two_tfs <- function()
{
  printf("--- test_1_footprint_two_tfs")

  loc <- tbl.fpAnnotated[grep("MA0679", tbl.fpAnnotated$motifName)[1], c("chr", "mfpStart", "mfpEnd")]
  target.tf <- c("SATB2", "CUX1")
  tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                             chromosome=loc$chr, startLoc=loc$mfpStart, endLoc=loc$mfpEnd,
                             transcriptionFactors=target.tf)

    # tbl.fpAnnotated identifies 3 motifs in this one footprint, and each of them is associated with SATB2
    # we therefore expect 3 rows, each describing the same footprint, and each with a different motif

  checkEquals(dim(tbl), c(3, 18))
  checkEquals(tbl$tfsMatched, rep("SATB2;CUX1", 3))
  checkEquals(tbl$motifName, c("MA0679.1", "MA0756.1", "MA0757.1"))

} # test_1_footprint_two_tfs
#------------------------------------------------------------------------------------------------------------------------
test_emptyResults <- function()
{
   printf("--- test_emptyResults")

      # first test out a chromosomal region 1 base long

   loc <- list(chr="chr1", start=3200, end=3200)
   target.tf <- "ELF3"   # mapped to the largest number of motifs in tbl.motifToMultipleGenes
   tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                          chromosome=loc$chr, startLoc=loc$start, endLoc=loc$end,
                          transcriptionFactors=target.tf)
   checkEquals(nrow(tbl), 0)

      # now repeat the SATB2 search (see above) with locs just 1 base too small, by incrementing start

   loc <- tbl.fpAnnotated[grep("MA0679", tbl.fpAnnotated$motifName)[1], c("chr", "motifStart", "motifEnd")]
   loc.orig <- loc
   loc$motifStart <- loc$motifStart + 1
   target.tf <- "SATB2"
   tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                          chromosome=loc$chr, startLoc=loc$motifStart, endLoc=loc$motifEnd,
                          transcriptionFactors=target.tf)
   checkEquals(nrow(tbl), 0)

      # now repeat the SATB2 search (see above) with locs just 1 base too small, by decrementing end

   loc <- loc.orig
   loc$motifEnd <- loc$motifEnd - 1
   target.tf <- "SATB2"
   tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                          chromosome=loc$chr, startLoc=loc$start, endLoc=loc$end,
                          transcriptionFactors=target.tf)
   checkEquals(nrow(tbl), 0)


} # test_emptyResults
#------------------------------------------------------------------------------------------------------------------------
test_360k_chr5_alzheimers_region <- function()
{
   printf("--- test_360k_chr5_alzheimers_region")

   loc <- list(chr="chr5", start=88640561, end=89001975)
   target.tfs <- c("POU5F2", "SATB2", "HLF", "SOX4", "CUX2")
   tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                              chromosome=loc$chr, startLoc=loc$start, endLoc=loc$end,
                              transcriptionFactors=target.tfs)
   checkEquals(dim(tbl), c(212, 18))
   checkEquals(sort(unique(tbl$tfsMatched)), c("HLF", "POU5F2", "POU5F2;SOX4", "SATB2;CUX2", "SOX4"))

} # test_360k_chr5_alzheimers_region
#------------------------------------------------------------------------------------------------------------------------
test_findSNPsInFootprints <- function()
{
   printf("--- test_findSNPsInFootprints")
   snp.chromosome <- "chr5"
   snp.loc.1 <- 88884587  # subset(tbl.snps, SNP %in% c("rs3814428", "rs770463", "rs446500"))$GRCh38.liftover
   snp.loc.3 <- 88883758  # from 1000 genomes, should retrieve 6-motif footprint

   tfs.all <- sort(unique(tbl.motifToMultipleGenes$tfs))
   tfs <- c("ZNF263", "MAZ", "KLF16")

       # search region
   region.chromosome <- "chr5"
   region.startLoc <- 88884570
   region.endLoc   <- 88884600   # chr5:88884570-88884600

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               region.chromosome, region.startLoc, region.endLoc,
                               snp.chromosome, snp.loc.1,
                               transcriptionFactors=tfs)

   checkEquals(dim(tbl), c(7, 9))
   checkEquals(sort(unique(tbl$tfsMatched)), c("KLF16", "MAZ", "ZNF263"))

      # now retrieve a snp in a footprint, but notice that it misses the reported motif sequence
   snp.loc.2 <- 88899133  # subset(tbl.snps, SNP %in% c("rs3814428", "rs770463", "rs446500"))$GRCh38.liftover
   region.startLoc <- 88899100
   region.endLoc   <- 88899300
   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               region.chromosome, region.startLoc, region.endLoc,   # region of interest
                               snp.chromosome, snp.loc.2,
                               transcriptionFactors="TCF7")
   checkEquals(dim(tbl), c(1,9))
   x <- as.list(tbl[1,])
   checkEquals(x$chr, "chr5")
   checkEquals(x$mfpStart, 88899109)
   checkEquals(x$mfpEnd, 88899140)
   checkEquals(x$motifStart, 88899109)
   checkEquals(x$motifEnd, 88899120)
   checkEquals(x$sequence, "AAAGTTCAAAGC")
   checkEquals(x$motifName, "MA0769.1")
   checkEquals(x$snp, 88899133)
   checkEquals(x$tfsMatched, "TCF7")

       # now try combining the two snp.locs successfully tested above
   tfs <- c("KLF16", "MAZ", "ZNF263", "TCF7")
   region.startLoc <- 88884570
   region.endLoc   <- 88899300

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               region.chromosome, region.startLoc, region.endLoc,   # region of interest
                               snp.chromosome, c(snp.loc.1, snp.loc.2),
                               transcriptionFactors=tfs)

   checkEquals(dim(tbl), c(8, 9))
   checkEquals(sort(unique(tbl$tfsMatched)), c("KLF16", "MAZ", "TCF7", "ZNF263"))
   checkEquals(sort(unique(tbl$snp)), c(88884587, 88899133))


} # test_findSNPsInFootprints
#------------------------------------------------------------------------------------------------------------------------
# snp rs10038371 is just upstream  chr5:88,903,407 of a footprint+motif at 88,903,378-403 (upstream from the perspective of the
# possibly regulated gene, MEF2C, which is encoded on the minus strand).
# here we look for that snp near that footprint, passing in a padding parameter which produces a liberal
# notion of overlapping
test_findSNPsNearFootprints <- function()
{
   printf("--- test_findSNPsNearFootprints")

   snp.chromosome <- "chr5"
   snp.loc <- 88903407
   target.region.chromosome <- "chr5"
   target.region.startLoc <- 88903000
   target.region.endLoc   <- 88903500

   tfs.all <- sort(unique(tbl.motifToMultipleGenes$tfs))
   #tfs <- c("KLF16", "MAZ", "ZNF263", "TCF7")

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               target.region.chromosome, target.region.startLoc, target.region.endLoc,
                               snp.chromosome, snp.loc,
                               transcriptionFactors=NA, # tfs.all,
                               padding=0)
   checkEquals(nrow(tbl), 0)

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               target.region.chromosome, target.region.startLoc, target.region.endLoc,
                               snp.chromosome, snp.loc,
                               transcriptionFactors=NA, # tfs.all,
                               padding=5)
   checkEquals(dim(tbl), c(1,9))

     # tfsMatched in the above table: ISL1;ISL2;LHX1;LHX2;LHX3;LHX4;LHX5;LHX6;LHX8;LHX9

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               target.region.chromosome, target.region.startLoc, target.region.endLoc,
                               snp.chromosome, snp.loc,
                               transcriptionFactors="ISL1",
                               padding=5)
   checkEquals(dim(tbl), c(1, 9))

     # with very large padding, our snp rs10038371 "overlaps" with the footprint at 88,902,507:88,902534
     # though of no likely biological interest, this helps to test the padding calculation
   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               target.region.chromosome, target.region.startLoc-1000, target.region.endLoc,
                               snp.chromosome, snp.loc,
                               transcriptionFactors=NA, #tfs.all,
                               padding=1000)

   checkEquals(dim(tbl), c(2, 9))

} # test_findSNPsNearFootprints
#------------------------------------------------------------------------------------------------------------------------
# a cryptic and transient bug, now no longer evident.
# originally, with all possible transcription factors considered, LHX1 was found (along with others in its gene family)
# for the region looked at here.
# but when LHX1 was considered alone, the results came back empty.
# this function ensures that they come back properly: one row describing an intersecting snp/fp on chr5 just upstream
# on the minus strand from mef2c
test_LXH1_matching <- function()
{
   printf("--- test_LXH1_matching")

   snp.chromosome <- "chr5"
   snp.loc <- 88903407
   target.region.chromosome <- "chr5"
   target.region.startLoc <- 88903000
   target.region.endLoc   <- 88903500

   tfs.all <- sort(unique(tbl.motifToMultipleGenes$tfs))

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               target.region.chromosome, target.region.startLoc, target.region.endLoc,
                               snp.chromosome, snp.loc,
                               transcriptionFactors=NA, # tfs.all,
                               padding=5)

   checkEquals(dim(tbl), c(1,9))
   target.tf <- "LHX1"
   checkTrue(grepl(target.tf, tbl[1, "tfsMatched"]))
   tokens <- strsplit(tbl[1, "tfsMatched"], ";")[[1]]
   checkTrue(match(target.tf, tokens) > 0)

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               target.region.chromosome, target.region.startLoc, target.region.endLoc,
                               snp.chromosome, snp.loc,
                               transcriptionFactors=target.tf,
                               padding=5)
   checkEquals(dim(tbl), c(1,9))

} # test_LHX1_matching
#------------------------------------------------------------------------------------------------------------------------
explore_abca7_snps <- function()
{
    printf("--- explore_abca7_snps")
    f <- system.file(package="PrivateCoryData", "data", "abca7-associated-snps-1kg.tsv")
    stopifnot(file.exists(f))
    tbl.snps <- read.table(f, header=FALSE, sep="\t", as.is=TRUE)
    colnames(tbl.snps) <- c("chrom", "loc")
    displayBedTable(igv, tbl.snps[, c(1,2,2)], "1kg snps")
    min.loc <- min(tbl.snps$loc)
    max.loc <- max(tbl.snps$loc)
    tbl.sub <- subset(tbl.fpAnnotated, chr=="chr19" & mfpStart >= (min.loc-500) & mfpEnd <= (max.loc + 1600))
    displayBedTable(igv, tbl.sub[, c("chr", "mfpStart", "mfpEnd")], "footprints")

   tbl <- findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                               "chr19", min.loc, max.loc,
                               "chr19", sort(unique(tbl.snps$loc)),
                               transcriptionFactors=NA,
                               padding=10)

} # explore_abca7_snps
#------------------------------------------------------------------------------------------------------------------------
test_snpFootDisplay <- function()
{
   printf("--- test_snpFootDisplay")
   if(!exists("igv"))
      igv <- igvR()

   checkTrue(connected(igv))

   #snpFootDisplay("chr5", start.loc, end.loc, genome="hg38")

   vcfDirectory <- system.file(package="igvR", "extdata", "vcf")
   vcfFile <- file.path(vcfDirectory, "chr22-sub.vcf.gz")
   tbifile <- file.path(vcfDirectory, "chr22-sub.vcf.gz.tbi")

   checkTrue(file.exists(vcfFile))
   checkTrue(file.exists(tbiFile))

   checkTrue(connected(igv))

   samples.01 <- c("HG00096", "HG00097", "HG00099")
   samples.02 <- c("HG00100", "HG00101")

   start.loc <- 49902169
   end.loc <- 49920831

   goto(igv, "chr22", start.loc, end.loc)
   displayVcfRegion(igv, "chr22", start.loc, end.loc, vcfFile, sampleIDs=samples.01)
   displayVcfRegion(igv, "chr22", start.loc, end.loc, vcfFile,  sampleIDs=samples.02)
   displayVcfRegion(igv, "chr22", start.loc, end.loc, vcfFile)


} # test_snpFootDisplay
#------------------------------------------------------------------------------------------------------------------------
