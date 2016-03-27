#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.onLoad <- function(libname, pkgname)
{
   printf("snpFoot package initializations:");
   printf("   loading tbl.fpAnnotated")
   data(tbl.fpAnnotated)  # may take 30 seconds or so: > 4M DNAse I footprints with TF motifs

   printf("   loading %s", "tbl.motifToMultipleGenes")
   data(tbl.motifToMultipleGenes)  # much quicker: maps motif identifiers to gene symbols (from Seth)
   printf("   loading tbl.gwasADsnpsInFp.05pval.igap2013")
   data(tbl.gwasADsnpsInFp.05pval.igap2013)

   printf("   loading tbl.gwas.level_1 AD, igap 2013")
   data(tbl.gwas.level_1)

   printf("   loading tbl.humangene3877")
   data(tbl.humangene3877)

} # .onLoad
#------------------------------------------------------------------------------------------------------------------------
