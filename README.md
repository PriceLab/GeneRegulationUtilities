# GeneRegulationUtilities
An R package with miscellaneous functions and public data for the exploration and analysis of gene regulation

_findTFsInFootprints_: for selected transcription factors, find  DNase I footprints in which they may bind

    library(GeneRegulationUtiltities)
    data(tbl.fpAnnotated)  # may take 30 seconds or so: > 4M DNAse I footprints with TF motifs
    data(tbl.motifToMultipleGenes  # much quicker: maps motif identifiers to gene symbols (from Seth)
    
    target.tfs <- c("POU5F2", "SATB2", "HLF", "SOX4", "CUX2")
    tbl.out <- findTFsInFootprints(tbl.fpAnnotated, tbl.motifToMultipleGenes, "chr5", 88640561, 89001975, target.tfs)
    dim(tbl.out)   # 212 17
    head(tbl.out)[, c(1:4, 17)]
            chr      start      end score tfsMatched
    3610411 chr5 88643612 88643621  53.5       SOX4
    3610412 chr5 88643613 88643621  50.5       SOX4
    3610451 chr5 88653378 88653392  52.6       SOX4
    3610452 chr5 88653378 88653392  61.2       SOX4
    3610456 chr5 88653569 88653578  52.7     POU5F2
    3610497 chr5 88656492 88656507  51.1     POU5F2

