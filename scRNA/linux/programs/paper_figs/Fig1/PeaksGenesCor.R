library(ggplot2)
library(viridis)
library(Hmisc)
library(plyr)
library(Seurat)
library(Signac)
library(foreach)
library(doParallel)
library(GenomicRanges)

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'MG', 'ODC', 'OPC', 'SUB')

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")

targetDir = "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/"
dir.create(targetDir, recursive = T, showWarnings = F)

getPeakRna = function(atacFilt, test, cluster, targetDir) {
  try({
    atacSub = subset(x = atacFilt, subset = Annotations == cluster)
    peakRna = LinkPeaks(
      object = atacSub,
      peak.assay = 'Combined_peaks',
      expression.assay = "RNA",
      peak.slot = "counts",
      expression.slot = "data",
      method = test,
      distance = 20000,
      min.cells = 10,
      pvalue_cutoff = 0.01,
      score_cutoff = 0.05,
      verbose = TRUE)
    p2g = data.frame(Links(peakRna))
    p2g$Cluster = cluster
    write.csv(p2g, file = paste0(targetDir, "peakGeneCor_", cluster, "_20K_2023-10-27.csv"), row.names = F)
    pritn(cluster)
  })
  return(p2g)
}


useCores<-2

cl <- makeCluster(useCores, type = "FORK")
registerDoParallel(cl)

p2gComb<-data.frame(foreach(i=clusters, .combine=rbind, .packages=c('Seurat', 'Signac') ) %dopar%{
  getPeakRna(atacFilt= atacFilt, test = "pearson", cluster = i, targetDir = "atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/")
})

parallel::stopCluster(cl = cl)


write.csv(p2gComb, file = paste0(targetDir, "peakGenes_All_20K_2023-10-27.csv"), row.names = F)

