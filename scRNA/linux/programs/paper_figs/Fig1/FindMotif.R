#BiocManager::install(version = "3.18")


library(Signac)
library(JASPAR2022)
library(Seurat)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
#library(patchwork)
library(JASPAR2020)

markers = read.csv("atacRna/peaksGenesCor/atacIntegrated_macs2_2_RNA/pearson/peakGeneCor_AllClusters_500K_2023-11-06_FiltClust.csv")

atacFilt = readRDS("atacIntegrated_macs2_2_RNA")
DefaultAssay(atacFilt) <- "Combined_peaks"
atacFilt[["PredictActivity"]] = NULL
gc()
# atacFilt <- FindTopFeatures(atacFilt, min.cutoff = 5)
# atacFilt <- RunTFIDF(atacFilt)
# 
# atacFilt<- RunSVD(atacFilt)


pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


atacFilt <- AddMotifs(
  object = atacFilt,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


all(unique(markers$peak) %in% rownames(atacFilt))

atacFilt <- RegionStats(atacFilt, genome = BSgenome.Mmusculus.UCSC.mm10)


curCluster = "CA1"
curMarkers = unique(markers$peak[markers$Cluster == curCluster])



enriched.motifs <- FindMotifs(
  object = atacFilt,
  features = curMarkers, assay= "Combined_peaks"
)

write.csv(enriched.motifs, "atacIntegrated_macs2_2_RNA_Motiffs_Jaspar2022.csv")



