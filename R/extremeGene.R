# extremeGene.R
#------------------------------------------------------------------------------------------------------------------------
library(biomaRt)
library(edgeR)
library(limma)
library(RUnit)
stringsAsFactors = FALSE;
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   dataDir <- "/mnt/data"
   dataDir <- "/Volumes/local/Cory/for_Paul/extreme"
   data <- loadData(dataDir)
   printf("--- testing threshold calculation")
   thresholds <- calculateHiLoThresholds(data$expr, "ENSG00000169245", 0.2)
   checkEquals(thresholds$hi, 35.6)
   checkEquals(thresholds$lo, 4)

} # runTess
#------------------------------------------------------------------------------------------------------------------------
loadData <- function(dataDir)
{

  f <- file.path(dataDir, "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates_Flowcell_1.csv")
  printf("loading metadata file '%s'", f)
  stopifnot(file.exists(f))
  tbl.meta  <- read.table(f, header = T, sep = ",")

  f <- file.path(dataDir, "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_JulyRerun.txt")
  stopifnot(file.exists(f))
  printf("loading expression file '%s'", f)
  mtx.expr <- as.matrix(read.table(f, header = T, sep = "", check.names = F, strip.white=TRUE))

  invisible(list(meta=tbl.meta, expr=mtx.expr))

} # loadData
#------------------------------------------------------------------------------------------------------------------------
calculateHiLoThresholds <- function(mtx.expr, gene.ensgID, threshold=0.20)
{
   stopifnot(gene.ensgID %in% rownames(mtx.expr))
   hi <- quantile(mtx.expr[gene.ensgID,], 1.0 - threshold)
   lo <- quantile(mtx.expr[gene.ensgID,], threshold)

   return(list(hi=as.numeric(hi), lo=as.numeric(lo)))

} # calculateHiLoThresholds
#------------------------------------------------------------------------------------------------------------------------
