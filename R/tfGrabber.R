library(snpFoot)
library(hash)
library(TReNA)
library(PrivateCoryData)
library(igvR)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv"))
    igv <- igvR()
#------------------------------------------------------------------------------------------------------------------------
TF_grabber<-function(genelist,trn.list,label=NULL,promoterDist=1000000)
{
   printf(" --- entering TF_grabber for %s", label)

 # input paras:
 # gene list: a vector of genes such as MEF2C, ABCA7, CLU. In the format of RefSeq gene symbol

  #genelist=c("MEF2C","ABCA7","CR1")
  #genelist=rownames(trn.rtrim_isb_AD)
  # requirements
  #require(hash)
  #require(GeneRegulationUtilities)
  #require(snpFoot)
  if (!exists("tbl.motifToMultipleGenes")){data(tbl.motifToMultipleGenes)}
  if (!exists("tbl.fpAnnotated")){data(tbl.fpAnnotated)}
  if (!exists("humangene")) {data(tbl.humangene3877)}

    #label of each TRN
  #browser()
  ntrn=length(trn.list)
  if (is.null(label)){label=paste("trn",1:ntrn,sep="")}

  # all TFs
  alltfs=vector()
  for (i in 1:ntrn){
    alltfs=c(alltfs,colnames(trn.list[[i]]))
    }
  alltfs=unique(alltfs)
  nalltfs=length(alltfs)

  # initialize TRN matrix
  singleTRN0=as.data.frame(matrix(0,nrow=length(trn.list),ncol=nalltfs))
  colnames(singleTRN0)=alltfs
  rownames(singleTRN0)=label

   ##################################################
   ########## main loop for each gene ###############
   ##################################################
   allTRN=list()
   wholebed=data.frame()
   for (gene in genelist){
      singleTRN=singleTRN0

     for (k in 1:length(trn.list)){# get single and summed TRN
        itrn=as.matrix(trn.list[[k]])
        if (!gene %in% rownames(itrn)){
           next;
           }
        else{
           subtrn=subset(itrn,rownames(itrn)==gene)
           singleTRN[k,colnames(subtrn)]=subtrn
           } # else
      } # for k

      # clean up single TRN

     #keeperCols <- which(colSums(abs(singleTRN)) > 0.01)
     #keeperRows <- which(rowSums(abs(singleTRN)) < 0.01)

         # remove all rows and columns with no values above threshold
     singleTRN=singleTRN[,which(colSums(abs(singleTRN)<0.01)!=ntrn)]
     singleTRN=singleTRN[which(rowSums(abs(singleTRN) <0.01)!=ncol(singleTRN)),]

# sumTRN
if (min(dim(singleTRN)) == 1){
  sumTRN=singleTRN
}else{
  sumTRN=as.data.frame(matrix(colSums(singleTRN),nrow=1))
  colnames(sumTRN)=colnames(singleTRN)
}


# get motif of each TF in the sumTRN
imotifs=subset(tbl.motifToMultipleGenes,tfs %in% colnames(sumTRN))
geneline=humangene[which(humangene$genename==gene),][1,]
ichr=as.character(geneline[,"chrom"])
startPosition=geneline[,"start"]-promoterDist
endPosition=geneline[,"start"]+promoterDist

tbl.thisgene <- subset(tbl.fpAnnotated,  chr==ichr & motifName %in% imotifs$motif & mfpStart >= startPosition & mfpEnd <= endPosition)
tbl.thisgene <- tbl.thisgene[, c("chr", "mfpStart", "mfpEnd", "motifName")]

if(nrow(tbl.thisgene) > 0 && nrow(imotifs) > 0){

# add regression coefficients
for (i in 1:nrow(tbl.thisgene)){
  thismotif=tbl.thisgene[i,"motifName"]
  #printf("motif (%d/%d): %s", i, nrow(tbl.thisgene), thismotif)
  itfs=subset(imotifs,motif==thismotif)
  info.string <- "<html>";
  gene.motif <- sprintf("%s:&nbsp;%s", gene, thismotif)
  info.string <- paste(info.string, gene.motif, sep="")
  for (j in 1:nrow(itfs)){
    #newname=paste(newname,paste(itfs[j,"tfs"],sumTRN[,itfs[j,"tfs"]],sep=":"),sep="-----")
    tf.info <- sprintf("<br>%s:&nbsp;%04.2f", itfs[j, "tfs"], sumTRN[,itfs[j,"tfs"]])
    #browser()
    info.string <- paste(info.string, tf.info, sep="")
    }
  info.string <- paste(info.string, "</html>", sep="")
  #print(info.string)
  tbl.thisgene[i,"name"]=info.string
  tbl.thisgene[i, "motifName"] <- thismotif
  } # for i

     #trn
   thistrn=list()
   thistrn$singleTRN=singleTRN
   thistrn$sumTRN=sumTRN
   allTRN[[gene]]=thistrn
   wholebed=rbind(wholebed,tbl.thisgene)
   } # if nrow

  }  # end of: for each gene in genelist

   #as.vector(as.list(h))
    TRN=list(TRN=allTRN,bed4igv=wholebed)
    return(TRN)

} # TF_grabber
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   # load the trns
   print(load("demoTRN.Rdata"))  # "trn.rtrim_isb_AD_cer" "trn.rtrim_isb_AD"
   dim(trn.rtrim_isb_AD_cer) # [1] 11822   500
   dim(trn.rtrim_isb_AD)     # [1] 8637  505
      # collect all TRNs into a list
   trn.list <- list(trn.rtrim_isb_AD_cer, trn.rtrim_isb_AD)
       # integrate TRN with footprint data: only 10 genes for demo
   #trn=TF_grabber(rownames(trn.rtrim_isb_AD)[2],trn.list,label=c("isb_AD_cer","isb_AD"),promoterDist = 1000000)
   trn=TF_grabber(rownames(trn.rtrim_isb_AD_cer)[1121],trn.list,label=c("isb_AD_cer","isb_AD"),promoterDist = 1000000)
      # display TRN in IGV
   connected(igv)
   displayBedTable(igv, trn$bed4igv,"test_TRN_track")
   #goto(igv, "chr1", 2986001, 248956422)  # FGR
   goto(igv, "chr5", 88544489, 89077857)   # MEF2C
   invisible(trn)

} # runTests
#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   dir <- "/Volumes/local/Cory/for_Hongdong"
   stopifnot(file.exists(dir))

   files <- grep("RData$", list.files(dir), value=TRUE)
   for(file in files[4:length(files)]){
      short.name <- sub(".RData", "", sub("^trn.rtrim_", "", file))
      full.path <- file.path(dir, file)
      stopifnot(file.exists(full.path))
      printf("--- loading %s", full.path)
          # learn the name of the trn-specific trimmed result variable, which is of the
          # form 'trn.rtrim.xxx', where xxx indicates its origin.
          # assign that specific variable name to the standard one, trn.rtrim, for use
          # below
      var.names <- load(full.path, envir=.GlobalEnv)
      var.name.trn.rtrim <- grep("^trn.rtrim", var.names, value=TRUE)
      eval(parse(text=sprintf("trn.rtrim <<- %s", var.name.trn.rtrim)))
      genes.of.interest <- intersect(rownames(trn.rtrim), c("MEF2C","ABCA7","CR1"))
      #genes.of.interest <- rownames(trn.rtrim)
      promoterDistance <- 10000
      printf("calling TF_grabber for '%s': %d genes, trn of dimension %d, %d, promoterDistance: %d",
             short.name, length(genes.of.interest), nrow(trn.rtrim), ncol(trn.rtrim), promoterDistance)
      time.info <- system.time(x <- TF_grabber(genes.of.interest, list(trn.rtrim), label=short.name, promoterDist=promoterDistance))
      output.filename <- sprintf("%s.%d-dist.%d-genes.results.RData", short.name, promoterDistance,
                                 length(genes.of.interest))
      printf("%10s: %10.2f seconds.  writing %d footprints on %d chromosomes to %s",
             short.name, time.info[["elapsed"]],
             nrow(x$bed4igv), length(unique(x$bed4igv$chr)), output.filename)
      save(x, file=output.filename)
      #browser()
      #zz <- 99
      }

} # run
#------------------------------------------------------------------------------------------------------------------------
# mariette.overlaps <- function()
# {
#    file <- "~/s/work/priceLab/cory/footprintFinderExperiments/mef2c-related-snps-from-mariette.tsv"
#    tbl.snps <- read.table(file, header=TRUE, as.is=TRUE, sep="\t")
#    snp.locs <- tbl.snps$GRCh38.liftover
#
#    print(range(snp.locs))
#    tbl.snpsFP <- findSNPsInFootprints(trn$bed4igv, tbl.motifToMultipleGenes,
#                                       "chr5", 88544489, 89077857,
#                                       "chr5", snp.locs, padding=1000)
#
# } # mariette.overlaps
#------------------------------------------------------------------------------------------------------------------------





