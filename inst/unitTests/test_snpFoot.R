library(RUnit)
library(snpFoot);
library(PrivateCoryData)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")
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
   test_tfGrabber()
   test_tfGrabber.new()
   test_run.tfGrabber.new()
   test_intersectingLocs()
   test_displayGWAS()
   test_intersectMarietSnpsWithGWAS()
   test_displayAllEncodeFootprints()
   test_runLift.hg19to38()
   test_selectBestAmongDuplicates()

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
# "ITGA2B" is a problematic gene in /Volumes/local/Cory/for_Hongdong/trn.rtrim_isb_all.RData
test_tfGrabber <- function()
{
   printf("--- test_tfGrabber")
   dir <- system.file(package="snpFoot", "extdata")
   checkTrue(file.exists(dir))
   file <- "demoTRN.RData"
   full.path <- file.path(dir, file)
   checkTrue(file.exists(full.path))
   var.names <- load(full.path, envir=.GlobalEnv)
   var.name.trn.rtrim <- grep("^trn.rtrim", var.names, value=TRUE)
   eval(parse(text=sprintf("trn.rtrim <<- %s", var.name.trn.rtrim[1])))
   genes.of.interest <- intersect(rownames(trn.rtrim), c("MEF2C","ABCA7","CR1"))
   promoterDistance <- 1000000

   time.info <- system.time(x <- tfGrabber(genes.of.interest, trn.rtrim,
                                           trnName="demo", promoterDist=promoterDistance))
   print(time.info)
   checkTrue(nrow(x$bed4igv) >= 10)
   checkEquals(length(x$TRN), length(genes.of.interest))

} # test_tfGrabber
#------------------------------------------------------------------------------------------------------------------------
test_tfGrabber.new <- function()
{
   printf("--- test_tfGrabber.new")

   dir <- system.file(package="snpFoot", "extdata")
   checkTrue(file.exists(dir))
   file <- "demoTRN.RData"
   full.path <- file.path(dir, file)
   checkTrue(file.exists(full.path))
   var.names <- load(full.path, envir=.GlobalEnv)
   var.name.trn.rtrim <- grep("^trn.rtrim", var.names, value=TRUE)
   eval(parse(text=sprintf("trn.rtrim <<- %s", var.name.trn.rtrim[1])))

   genes.of.interest <- intersect(rownames(trn.rtrim), c("MEF2C","ABCA7","CR1"))
   promoterDistance <- 1000000

   x <- tfGrabber.new("MEF2C", trn.rtrim, trnName="demo", promoterDist=10)
   checkEquals(dim(x), c(0,0))

   x <- tfGrabber.new("MEF2C", trn.rtrim, trnName="demo", promoterDist=10000)
   checkEquals(dim(x), c(4, 23))

   x <- tfGrabber.new("MEF2C", trn.rtrim, trnName="demo", promoterDist=1000000)
   checkEquals(dim(x), c(605, 23))

   x <- tfGrabber.new("CR1", trn.rtrim, trnName="demo", promoterDist=1000000)
   checkEquals(dim(x), c(0,0))

   xx <- lapply(c("MEF2C", "ABCA7"), function(gene) tfGrabber.new(gene, trn.rtrim, trnName="demo", promoterDist=1e6))
   xxx <- do.call("rbind", xx)

   genes.44 <-  c("ABCA7", "APOE", "BIN1", "CASS4", "CD2AP", "CELF1", "CLU",
                  "CR1", "EPHA1", "FERMT2", "HLA", "INPP5D", "MEF2C", "MS4A6A",
                  "NME8", "PICALM", "PTK2B", "SLC24A4/RIN3", "SORL1", "ZCWWPW1",
                  "ABI3", "IRF8", "MS4A14", "MS4A4A", "MS4A7", "PLCG1", "PLCG2",
                  "TGFB1", "TREM2", "TTC3", "UNC5C", "PILRB", "SLC24A4",
                  "FASLG", "TNFSF18", "PIEZO2", "NAPG", "APCDD1", "DNMBP",
                  "APP", "BACE1", "BACE2", "PSEN1", "PSEN2")
     # 5 seconds elapsed
   time.info <- system.time(fp.44.list <- lapply(genes.44, function(gene) tfGrabber.new(gene, trn.rtrim, trnName="demo", promoterDist=1e6)))
   lapply(fp.44.list, length)
   tbl.fp.44 <- do.call("rbind", fp.44.list)
   dim(tbl.fp.44) # [1] 8470   23
   table(tbl.fp.44$gene)

     #  ABCA7    ABI3  APCDD1    APOE   BACE1   BACE2    BIN1   CASS4   CD2AP   CELF1   DNMBP   EPHA1
     #    233     183     432     676    1117      57      18      62     167     177      57      30
     # FERMT2  INPP5D   MEF2C  MS4A14  MS4A4A  MS4A6A    NAPG  PIEZO2   PLCG1   PLCG2   PSEN1   PSEN2
     #    411     185     605     122      16      49     445     181     121     158     453     177
     #  PTK2B SLC24A4   SORL1   TGFB1   TREM2    TTC3   UNC5C
     #    300      84     218    1435      46     148     107

} # test_tfGrabber.new
#------------------------------------------------------------------------------------------------------------------------
test_run.tfGrabber.new <- function()
{
   printf("--- test_run.tfGrabber.new")
   load(system.file(package="PrivateCoryData", "extdata/trn10.genes44/genes.44.RData"))
   checkTrue(exists("genes.44"))
   checkEquals(length(genes.44), 44)

   dir <- system.file(package="snpFoot", "extdata")
   checkTrue(file.exists(dir))

   file <- "demoTRN.RData"
   full.path <- file.path(dir, file)
   checkTrue(file.exists(full.path))
   var.names <- load(full.path, envir=.GlobalEnv)
   var.name.trn.rtrim <- grep("^trn.rtrim", var.names, value=TRUE)
   eval(parse(text=sprintf("trn.rtrim <<- %s", var.name.trn.rtrim[1])))
   checkEquals(dim(trn.rtrim), c(11822, 500))

   genes.missing <- setdiff(genes.44, rownames(trn.rtrim))
   genes.missing
      # [1] "APP"     "CLU"     "CR1"     "FASLG"   "HLA"     "IRF8"    "MS4A7"   "NME8"    "PICALM"  "PILRB"   "TNFSF18" "ZCWWPW1"
   intersect(genes.missing, colnames(trn.rtrim)) # [1] "IRF8"
   genes.of.interest <- intersect(rownames(trn.rtrim), genes.44)   # 32

   sqlite.fileName <- file.path(getwd(), "grabber.sqlite")
   if(file.exists(sqlite.fileName))
       unlink(sqlite.fileName)

   checkTrue(!file.exists(sqlite.fileName))
   sql.tableName <- "tfFootprints"
   tbl.fpDemo <- run.tfGrabber(genes.of.interest, trn.rtrim, "trnDemo", sqlite.fileName, sql.tableName)
   db <- dbConnect(dbDriver("SQLite"), sqlite.fileName)
   checkTrue(file.exists(sqlite.fileName))

   result <- dbGetQuery(db, sprintf("select count(*) from %s", sql.tableName))
   rowCount <- result[1,1]
   checkTrue(rowCount > 8000)
   checkTrue(rowCount < 9000)

   tbl <- dbGetQuery(db, sprintf("select * from %s limit 3", sql.tableName))
   checkEquals(dim(tbl), c(3, 24))
   checkEquals(sort(colnames(tbl)), sort(c("TF", "beta", "chr", "end", "experimentCount", "fpEnd", "fpScore",
                                           "fpStart", "gene", "id", "motifEnd", "motifFpOverlap", "motifName",
                                           "motifScore", "motifStart", "motifStrand", "name", "promoterDist",
                                           "pvalue", "qvalue", "row_names", "sequence", "start", "trnName")))

} # test_run.tfGrabber.new
#------------------------------------------------------------------------------------------------------------------------
test_intersectingLocs <- function()
{
   printf("--- test_intersectingLocs")

   chrom <- "chr1"
   tbl.bed <- data.frame(chrom=chrom, start=1, end=10, stringsAsFactors=FALSE)
   loc.good <- 6
   loc.bad  <- 16
   checkEquals(intersectingLocs(tbl.bed, chrom, loc.good), loc.good)
   checkEquals(intersectingLocs(tbl.bed, chrom, loc.bad), numeric(0))
   checkEquals(intersectingLocs(tbl.bed, chrom, loc.bad, padding=10), loc.bad)

} # test_intersectingLocs
#------------------------------------------------------------------------------------------------------------------------
test_displayGWAS <- function()
{
   printf("--- test_displayGWAS")

   checkTrue(exists("tbl.gwas"))
   if(!exists("igv"))
       return(FALSE)

    displayGWASTable(igv, tbl.gwas, "ADgwas2013")

} # test_displayGWAS
#------------------------------------------------------------------------------------------------------------------------
test_intersectMarietSnpsWithGWAS <- function()
{
   printf("--- test_intersectMarietSnpsWithGWAS")

   if(!exists("igv"))
       return(FALSE)

     # first, display all mariette snps

   file <- system.file(package="PrivateCoryData", "extdata", "mef2c-related-snps-from-mariette.tsv")
   checkTrue(file.exists(file))
   tbl.mariette <- read.table(file, sep="\t", header=TRUE, as.is=TRUE)
   mariette.snps <- tbl.mariette$BP
   displaySnps(igv, "chr5", mariette.snps, "mariette.snps")
   goto(igv, "chr5", 87850803, 88515362)

     # now load all the footprints from 10trns, 44 genes
   file <- system.file(package="PrivateCoryData", "extdata", "trn10.genes44", "44genes-10trns-bedTable.RData")
   checkTrue(file.exists(file))
   load(file)   # it's called tbl.out, with almost 18k rows
   displayBedTable(igv, tbl.out, "44genes.10trns")  # igv will ask if you want to create an index.  say yes
   squishTrack(igv, "44genes.10trns.bed")

      # create a track of intersection of gwas snps with 44genes.10trns footprints
   gwas.chr5.snps <- subset(tbl.gwas, CHR=="chr5")$BP
   gwas.chr5.snps.in.fp <- intersectingLocs(tbl.out, "chr5", gwas.chr5.snps, padding=10)  # just two
   displaySnps(igv, "chr5", gwas.chr5.snps.in.fp, "gwas.sps.in.fp")

      # now do the same with mariette snps
   mariette.snps.in.or.near.fp <- intersectingLocs(tbl.out, "chr5", mariette.snps, padding=50)  # also just two
   displaySnps(igv, "chr5", mariette.snps.in.or.near.fp, "mariette.sps.in.fp")

} # test_intersectMarietSnpsWithGWAS
#------------------------------------------------------------------------------------------------------------------------
test_displayAllEncodeFootprints <- function()
{
   printf("--- test_displayAllEncodeFootprints")

   if(!exists("igv"))
      return()

   checkTrue(exists("tbl.fpAnnotated"))
   tbl.chr5 <- subset(tbl.fpAnnotated, chr=="chr5")
   displayBedTable(igv, tbl.chr5[, c("chr", "mfpStart", "mfpEnd", "motifName")], "from.encode")

} # test_displayAllEncodeFootprints
#------------------------------------------------------------------------------------------------------------------------
test_runLift.hg19to38 <- function()
{
    printf("--- test_runLift.hg19to38")

    small.bed.file <- system.file(package="PrivateCoryData", "extdata/elizabethBlue/chr1_CU0039F.shared.bed")

    expected.output.file <- "./chr1_CU0039F.shared.hg38.bed"
    if(file.exists("expected.output.file"))
       unlink(expected.output.file)

    runLift.hg19to38(small.bed.file)
    checkTrue(file.exists(expected.output.file))
    tbl <- read.table(expected.output.file, sep="\t", as.is=TRUE)
    checkEquals(dim(tbl), c(30, 4))
    checkEquals(tbl[1,1], "chr1")
    checkEquals(tbl[1,2], 172863016)
    checkEquals(tbl[1,3], 172863016)
    checkEquals(tbl[1,4], "1:172832156")

    if(file.exists("expected.output.file"))
       unlink(expected.output.file)

} # test_runLift.hg19to38
#------------------------------------------------------------------------------------------------------------------------
test_selectBestAmongDuplicates <- function()
{
   printf("--- test_selectBestAmongDuplicates")
   file.path <- system.file(package="snpFoot", "extdata", "footprints.tbl.withDups.RData")
   load(file.path)
   tbl.best <- selectBestAmongDuplicates(tbl.fpoi.snp)
   checkEquals(nrow(tbl.best), 6)

} # test_selectBestAmongDuplicates
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()) runTests()
