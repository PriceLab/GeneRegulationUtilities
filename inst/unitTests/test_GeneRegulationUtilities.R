library(RUnit)
library(GeneRegulationUtilities);
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.fpAnnotated")){
   printf("loading tbl.fpAnnotated from GeneRegulationUtilities package...")
   data(tbl.fpAnnotated)
   }

if(!exists("tbl.motifToMultipleGenes")){
   printf("loading tbl.motifToMultipleGenes from GeneRegulationUtilities package...")
   data(tbl.motifToMultipleGenes)
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_1_footprint_SAT2B_tf()
   test_1_footprint_two_tfs()
   test_emptyResults()
   test_360k_chr5_alzheimers_region()

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

  loc <- tbl.fpAnnotated[grep("MA0679", tbl.fpAnnotated$info)[1], c("chr", "start", "end")]
  target.tf <- "SATB2"
  tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                         chromosome=loc$chr, startLoc=loc$start, endLoc=loc$end,
                         transcriptionFactors=target.tf)

    # tbl.fpAnnotated identifies 3 motifs in this one footprint, and each of them is associated with SATB2
    # we therefore expect 3 rows, each describing the same footprint, and each with a different motif

  checkEquals(dim(tbl), c(3, 17))
  checkEquals(tbl$tfsMatched, rep("SATB2", 3))
  checkEquals(tbl$name, c("MA0679.1", "MA0756.1", "MA0757.1"))

} # test_1_footprint_SAT2B_tf
#------------------------------------------------------------------------------------------------------------------------
test_1_footprint_two_tfs <- function()
{
  printf("--- test_1_footprint_two_tfs")

  loc <- tbl.fpAnnotated[grep("MA0679", tbl.fpAnnotated$info)[1], c("chr", "start", "end")]
  target.tf <- c("SATB2", "CUX1")
  tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                         chromosome=loc$chr, startLoc=loc$start, endLoc=loc$end,
                         transcriptionFactors=target.tf)

    # tbl.fpAnnotated identifies 3 motifs in this one footprint, and each of them is associated with SATB2
    # we therefore expect 3 rows, each describing the same footprint, and each with a different motif

  checkEquals(dim(tbl), c(3, 17))
  checkEquals(tbl$tfsMatched, rep("SATB2;CUX1", 3))
  checkEquals(tbl$name, c("MA0679.1", "MA0756.1", "MA0757.1"))

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

   loc <- tbl.fpAnnotated[grep("MA0679", tbl.fpAnnotated$info)[1], c("chr", "start", "end")]
   loc.orig <- loc
   loc$start <- loc$start + 1
   target.tf <- "SATB2"   # mapped to the largest number of motifs in tbl.motifToMultipleGenes
   tbl <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                          chromosome=loc$chr, startLoc=loc$start, endLoc=loc$end,
                          transcriptionFactors=target.tf)
   checkEquals(nrow(tbl), 0)

      # now repeat the SATB2 search (see above) with locs just 1 base too small, by decrementing end

   loc <- loc.orig
   loc$end <- loc$end - 1
   target.tf <- "SATB2"   # mapped to the largest number of motifs in tbl.motifToMultipleGenes
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
   checkEquals(dim(tbl), c(212, 17))
   checkEquals(sort(unique(tbl$tfsMatched)), c("HLF", "POU5F2", "POU5F2;SOX4", "SATB2;CUX2", "SOX4"))

} # test_360k_chr5_alzheimers_region
#------------------------------------------------------------------------------------------------------------------------
