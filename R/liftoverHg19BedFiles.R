library(GenomicRanges)
library(rtracklayer)
library(PrivateCoryData)
#------------------------------------------------------------------------------------------------------------------------
chain.19to38 <- import.chain(system.file(package="PrivateCoryData", "data", "hg19ToHg38.over.chain"))
chain.38to19 <- import.chain(system.file(package="PrivateCoryData", "data", "hg38ToHg19.over.chain"))
#------------------------------------------------------------------------------------------------------------------------
liftoverBedFile.19.38 <- function(filename)
{
   tbl <- read.table(filename, sep="\t", as.is=TRUE)

     # some minimal sanity checks:  chrX or chromX or  X values in 1st column

   chromValues <- sort(unique(sub("^chr", "", tbl[,1])))
   chromValues <- sub("^chrom", "", chromValues)
   chromValues <- toupper(chromValues)
   expected.chrom.values <- all(chromValues %in% as.character(c(1:22, "X", "Y", "M")))
   if(!expected.chrom.values)
      stop(sprintf("unexpected chromosome names in input bed file: %s", paste(chromValues, collapse=",")))

     # all numeric in second and third columsn

   stopifnot(all(is.numeric(tbl[,2])))
   stopifnot(all(is.numeric(tbl[,3])))

   if(ncol(tbl) == 3)
       colnames(tbl) <- c("chrom", "start", "end")

   if(ncol(tbl) == 4){
      colnames(tbl) <- c("chrom", "start", "end", "name")
      }

   gr <- GRanges(Rle(tbl$chrom), IRanges(tbl$start, tbl$end))
   seqlevelsStyle(gr) <- "UCSC"
   seqinfo(gr) <- SeqinfoForUCSCGenome("hg19")[seqlevels(gr)]
   gr.38 <- IRanges::unlist(liftOver(gr, chain.19to38))
   seqinfo(gr.38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr)]
   tbl.out <- data.frame(chrom=as.character(seqnames(gr.38)), start=start(gr.38), end=end(gr.38))
   if("name" %in% colnames(tbl))
      tbl.out <- cbind(tbl.out, name=tbl$name)

   filename.base <- basename(filename)
   filename.out <- sub(".bed$", ".hg38.bed", filename.base)

   printf("writing %d rows to %s", nrow(tbl.out), file.path(getwd(), filename.out))

   write.table(tbl.out, file=filename.out, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

} # liftoverBedFile.19.38
#------------------------------------------------------------------------------------------------------------------------
runLift.hg19to38 <- function(hg19.bed.filenames)
{
  for(filename in hg19.bed.filenames){
     printf("can we read '%s'? %s", filename, file.exists(filename))
     stopifnot(file.exists(filename))
     }

   for(file in hg19.bed.filenames)
     liftoverBedFile.19.38(file)

} # runLift.hg19to38
#------------------------------------------------------------------------------------------------------------------------
