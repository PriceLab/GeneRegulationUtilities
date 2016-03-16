#------------------------------------------------------------------------------------------------------------------------
tfGrabber1<-function(genelist,inputTRN,label=NULL,promoterDist=1000000, threshold=0.01)
{
   # input paras:
   # gene list: a vector of genes such as MEF2C, ABCA7, CLU. In the format of RefSeq gene symbol

   #genelist=c("MEF2C","ABCA7","CR1")
   #genelist=rownames(trn.rtrim_isb_AD)


   # requirements
   #require(hash)
   #require(GeneRegulationUtilities)

  if (!exists("tbl.motifToMultipleGenes")){data(tbl.motifToMultipleGenes)}
  if (!exists("tbl.fpAnnotated")){data(tbl.fpAnnotated)}
  if (!exists("humangene")) {data(tbl.humangene3877)}

  #label of each TRN
  if (is.null(label)){label="trn"}

  # all TFs
  alltfs=colnames(inputTRN)
  allgenes=rownames(inputTRN)
  nalltfs=length(alltfs)
  genelist=intersect(genelist,allgenes)
  if (length(genelist) == 0){warning("Input genes not in TRN")}

  allTRN=list()
  wholebed=data.frame()

  count <- 0
  for (gene in genelist){     # get single and summed TRN
    count <- count + 1
    if (!gene %in% rownames(inputTRN)) next;
    singleTRN=inputTRN[gene,]
    potential.regulator.count <- length(singleTRN)
    keepCols=which(abs(singleTRN)>threshold)     # keep only TFs with absolute coefficients above threshold
    singleTRN=singleTRN[keepCols]
    printf("--- %s %d) %8s regulators above threshold %d/%d", label, count, gene,  length(singleTRN), potential.regulator.count)
    # sumTRN == singleTRN in thse case of input is only one TRN
    sumTRN=singleTRN

    # get motif of each TF in the sumTRN
    imotifs=subset(tbl.motifToMultipleGenes,tfs %in% names(sumTRN))
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
         itfs=subset(imotifs,motif==thismotif)
         info.string <- "<html>";
         gene.motif <- sprintf("%s:&nbsp;%s", gene, thismotif)
         info.string <- paste(info.string, gene.motif, sep="")
         for (j in 1:nrow(itfs)){
            tf.info <- sprintf("<br>%s:&nbsp;%04.2f", itfs[j, "tfs"], sumTRN[itfs[j,"tfs"]])
            info.string <- paste(info.string, tf.info, sep="")
            }
         info.string <- paste(info.string, "</html>", sep="")
         tbl.thisgene[i,"name"]=info.string
         tbl.thisgene[i, "motifName"] <- thismotif
         } # for i
      thistrn=list()
      thistrn$singleTRN=singleTRN
      thistrn$sumTRN=sumTRN
      allTRN[[gene]]=thistrn
      wholebed=rbind(wholebed,tbl.thisgene)
    } # if nrow
  }  # for gene


  # return
  TRN=list(TRN=allTRN,bed4igv=wholebed)
  return(TRN)

} # tfGrabber1
#------------------------------------------------------------------------------------------------------------------------
