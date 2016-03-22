library(rtracklayer)
library(PrivateCoryData)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
chain.19to38 <- import.chain(system.file(package="PrivateCoryData", "data", "hg19ToHg38.over.chain"))
chain.38to19 <- import.chain(system.file(package="PrivateCoryData", "data", "hg38ToHg19.over.chain"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_liftoverBedFile.19.38()

} # runTest
#------------------------------------------------------------------------------------------------------------------------
liftoverBedFile.19.38 <- function(filename)
{
   tbl <- read.table(filename, sep="\t", header=TRUE, as.is=TRUE)
   gr <- GRanges(Rle(tbl$chrom), IRanges(tbl$start, tbl$end))
   seqlevelsStyle(gr) <- "UCSC"
   seqinfo(gr) <- SeqinfoForUCSCGenome("hg19")[seqlevels(gr)]
   gr.38 <- unlist(liftOver(gr, chain.19to38))
   seqinfo(gr.38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr)]
   tbl.out <- as.data.frame(gr.38)
   filename.out <- sub(".bed$", ".hg38.bed", filename)
   printf("writing %d rows to %s", nrow(tbl.out), filename.out)
   write.table(tbl.out[, c(1,2,3)], file=filename.out, col.names=FALSE, row.names=FALSE, quote=FALSE)

} # liftoverBedFile.19.38
#------------------------------------------------------------------------------------------------------------------------
test_liftoverBedFile.19.38 <- function ()
{
    print("--- test_liftoverBedFile.19.38")
    file <- "chr18_CU0070F.shared.bed"
    liftoverBedFile.19.38(file)

} # test_liftoverBedFile.19.38
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
    hg19.bed.files <- grep("shared.bed$", list.files("."), v=TRUE)
    for(file in hg19.bed.files)
      liftoverBedFile.19.38(file)

    hg38.bed.files <- grep("38.bed$", list.files("."), v=TRUE)
    if(!exists("igv"))
        igv < igvR()

    if(connected(igv))
       for(file in hg38.bed.files){
         loadFile(igv, file.path(getwd(), file))
         }

} # run
#------------------------------------------------------------------------------------------------------------------------

