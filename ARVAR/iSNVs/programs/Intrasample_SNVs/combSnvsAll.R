# Ampseq 

combDat = function(path) {
  combDat = data.frame(matrix(nrow = 0, ncol =0))
  samples = list.files(path)
  for (curSample in samples) {
    try({
      curDf = paste0(path, "/", curSample, "/filtered.csv")
      #print(curDf)
      if (!file.exists(curDf)) {
        print(list.files(paste0(path, "/", curSample)))
      }
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

ampseq = combDat(path = "IntraSnv_ampseq_all")

# save ampseq

write.csv(ampseq, "IntraSnv_results/ampseq_comb.csv", row.names = F)


metaseq = combDat(path = "IntraSnv_metaseq_all")

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
excludeSamples = unique(remCont(contamList, metaseq))
libsDf = metaseq[!metaseq$OrigName%in%excludeSamples,]

# save metaseq

write.csv(libsDf, "IntraSnv_results/metaseq_comb.csv", row.names = F)