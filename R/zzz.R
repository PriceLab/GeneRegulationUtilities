#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
.onLoad <- function(libname, pkgname)
{
   printf("snpFoot package initializations:");

   printf("   loading tbl.fpAnnotated")
       # may take 30 seconds or so: > 4M DNAse I footprints with TF motifs
   load(system.file(package="PrivateCoryData", "data", "tbl.fpAnnotated.RData"), envir=.GlobalEnv)

   printf("   loading %s", "tbl.motifToMultipleGenes")
   load(system.file(package="PrivateCoryData", "data", "tbl.motifToMultipleGenes.RData"), envir=.GlobalEnv)

   printf("   loading tbl.gwasADsnpsInFp.05pval.igap2013")
   load(system.file(package="PrivateCoryData", "data", "tbl.gwasADsnpsInFp.05pval.igap2013.RData"), envir=.GlobalEnv)

   printf("   loading tbl.gwas.level_1 AD, igap 2013")
   load(system.file(package="PrivateCoryData", "data", "tbl.gwas.level_1.RData"), envir=.GlobalEnv)

   printf("   loading tbl.humangene3877")
   load(system.file(package="PrivateCoryData", "data", "tbl.humangene3877.RData"), envir=.GlobalEnv)

} # .onLoad
#------------------------------------------------------------------------------------------------------------------------
