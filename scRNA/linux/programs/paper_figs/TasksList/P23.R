library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

rnaDat<-readRDS('integ_OPC_ODC')


renameClusters = function(df, clustCol) {
  newClust = character()
  for ( i in 1:nrow(df)) {
    if (df[i, clustCol] == 1) {
      curClust = "ODC"
    } else if (df[i, clustCol] == 2) {
      curClust = "OPC"
    } else if ((df[i, clustCol] == 3) | (df[i, clustCol] == 4)) {
      curClust = "Intermideate"
    }
    newClust = c(newClust, curClust)
  }
  df$newMonocClust = newClust
  return(df)
}

rnaDat@meta.data = renameClusters(df = rnaDat@meta.data, clustCol = "MonocClust")

opcOdcCells = rownames(rnaDat@meta.data)

addEnrich<-function(x){
  load(x)
  AUCmat <- AUCell::getAUC(cells_AUC)
  AucSub = AUCmat[, c(opcOdcCells)]
  rnaDat[['AUC']] <- CreateAssayObject(data = AucSub)
  DefaultAssay(rnaDat) <- 'AUC'
  rnaDat <- ScaleData(rnaDat, assay = 'AUC', features = rownames(AUCmat))
  return(rnaDat)
}

# add AUcell data
rnaDat<-addEnrich(x='./cellsAUC_keggClustProf.RData')
DefaultAssay(rnaDat)

colnames(rnaDat@meta.data)

Idents(rnaDat) = rnaDat$group

odcDat = subset(rnaDat, subset = newMonocClust == "ODC")

allMarkers <- FindMarkers(odcDat , only.pos = T, min.pct = 0.1, logfc.threshold = 0, test.use = "wilcox", ident.1 = "Stress", ident.2 = "Control")
allMarkers$Description = rownames(allMarkers)

fer = allMarkers[allMarkers$Description == "Ferroptosis",]