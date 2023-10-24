

# Ampseq 

combDat = function(path) {
  combDat = data.frame(matrix(nrow = 0, ncol =0))
  samples = list.files(path)
  for (curSample in samples) {
    try({
      curDf = paste0(path, "/", curSample, "/filtered.csv")
      df = read.csv(curDf)
      exactSample = gsub("_", "-", curSample)
      exactSample = sub(".*EHC", "EHC", exactSample)
      df$OrigName = curSample
      df$ExactSample = exactSample
      sampNamelist = strsplit(exactSample, "-")
      sampleName =  sampNamelist[[1]][1:3]
      sampleName = paste(sampleName, collapse = "-")
      df$Sample = sampleName
      df$Samp_Pos_Ref_Alt = paste(df$Sample, df$POSITION, df$REF.NT, df$VAR.NT, sep = "__")
      df$Path = curDf
      
      combDat = rbind(combDat, df)
    })
    
  }
  return(combDat)
}

ampseq = combDat(path = "IntraSnv_ampseq_overlap")

# save ampseq

write.csv(ampseq, "IntraSnv_ampseq_overlap/ampseq_comb.csv", row.names = F)

length(unique(ampseq$Sample)) == length(unique(ampseq$OrigName))

# metaseq 
metaseq = combDat(path = "IntraSnv_metaseq_overlap")
length(unique(metaseq$Sample)) == length(unique(metaseq$OrigName))
metaNames = unique(metaseq[, c("Sample", "OrigName")])
freqDf = data.frame(table((metaNames$Sample)))
selDf = metaseq[grepl("EHC-C19-1636Z", metaseq$Path),]
unique(selDf$REF.GENOME)
selSample = selDf[selDf$REF.GENOME == "p21175-s016_EHC-C19-1636Z_S16_L001",]
selStats = metaStats[grepl("EHC-C19-1636Z", metaStats$Sample_id),]
metaseqFilt = metaseq[!metaseq$OrigName == "p21175-s016_EHC-C19-1636Z_S16_L001",]

# remove contamination
remCont = function(contamList, libsDf) {
  exclSamples = character()
  contamList = gsub("_", ".", contamList)
  samplesList = libsDf$OrigName
  for (i in contamList) {
    curSamples = samplesList[grepl(i, samplesList)]
    exclSamples = c(exclSamples, curSamples)
  }
  return(exclSamples)
}

contamList = read.table("metaseqContam.list")
contamList = contamList$V1
excludeSamples = remCont(contamList, metaseqFilt)
libsDf = metaseqFilt[!metaseqFilt$OrigName%in%excludeSamples,]
write.csv(libsDf, "IntraSnv_metaseq_overlap/metaseq_comb.csv", row.names = F)
