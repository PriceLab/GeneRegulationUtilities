# GeneRegulationUtilities
An R package with miscellaneous functions and public data for the exploration and analysis of gene regulation

_findTFsInFootprints_: for selected transcription factors, find  DNase I footprints in which they may bind

    library(GeneRegulationUtilities)
    data(tbl.fpAnnotated)  # may take 30 seconds or so: > 4M DNAse I footprints with TF motifs
    data(tbl.motifToMultipleGenes)  # much quicker: maps motif identifiers to gene symbols (from Seth)

    target.tfs <- c("POU5F2", "SATB2", "HLF", "SOX4", "CUX2")
    tbl.out <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes, "chr5", 88640561, 89001975, target.tfs)
    dim(tbl.out)   # 212 18
    head(tbl.out)[, c(1:4, 18)]
            chr      start      end score tfsMatched
    3610411 chr5 88643612 88643621  53.5       SOX4
    3610412 chr5 88643613 88643621  50.5       SOX4
    3610451 chr5 88653378 88653392  52.6       SOX4
    3610452 chr5 88653378 88653392  61.2       SOX4
    3610456 chr5 88653569 88653578  52.7     POU5F2
    3610497 chr5 88656492 88656507  51.1     POU5F2


_findSNPsInFootprints_: find the footprints which overlap with the specified SNPS

    chrom <- "chr19"
    min.loc <- 1004844
    max.loc <- 1100074
    candidate.snps <- c(1007158, 1095102, 1095122, 1095630)
    findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                         chrom, min.loc, max.loc,
                         chrom, candidate.snps,
                         transcriptionFactors=NA,
                         padding=0)

      chr mfpStart  mfpEnd motifStart motifEnd        sequence     name     snp	tfsMatched
    chr19  1095102 1095131    1095117  1095131 GAGGCCCAGAGGTCG MA0728.1 1095122	NR2F6
    chr19  1095630 1095657    1095630  1095643  ACTTCGCCCCCGCC MA0162.2 1095630 EGR1;EGR2;...

Add padding to the search so that snps _near_ footprints are also retured:

   findSNPsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes,
                        chrom, min.loc, max.loc,
                        chrom, candidate.snps,
                        transcriptionFactors=NA,
                        padding=10)

      chr mfpStart  mfpEnd motifStart motifEnd        sequence     name     snp	tfsMatched
    chr19  1095062 1095094    1095062  1095076 CCCCCTCCCTTCCCC MA0516.1 1095102 EGR1;EGR2;...
    chr19  1095102 1095131    1095117  1095131 GAGGCCCAGAGGTCG MA0728.1 1095122 NR2F6
    chr19  1095630 1095657    1095630  1095643  ACTTCGCCCCCGCC MA0162.2 1095630 EGR1;EGR2;...

Create an igv track from data/tbl.gwasADsnpsInFp.05pval.igap2013.RData

````
data(tbl.gwasADsnpsInFp.05pval.igap2013)
displayBedTable(igv, tbl.gwasADsnpsInFp[, c("chr", "snp", "snp")], "gwasFPsnps")
````
