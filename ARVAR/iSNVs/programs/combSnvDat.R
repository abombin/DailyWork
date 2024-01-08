targetSamples = read.csv("sorting_tables/SeqSamples_Eval_Depth10.csv")
targetPattern = gsub("_", ".", targetSamples$Sample_id)
targetPattern = gsub("-", ".", targetPattern)

# Ampseq 

combDat = function(path) {
  combDat = data.frame(matrix(nrow = 0, ncol =0))
  samples = list.files(path)
  for (curSample in samples) {
    try({
      curDf = paste0(path, "/", curSample, "/filtered.csv")
      df = read.csv(curDf)
      exactSample = gsub("_", "-", curSample)
      eaxactSample = sub(".*EHC", "EHC", exactSample)
      df$OrigName = curSample
      df$ExactSample = exactSample
      sampNamelist = strsplit(exactSample, "-")
      sampleName =  sampNamelist[[1]][1:3]
      sampleName = paste(sampleName, collapse = "-")
      df$Sample = sampleName
      df$Samp_Pos_Ref_Alt = paste(df$Sample, df$POSITION, df$REF.NT, df$VAR.NT, sep = "__")
      
      combDat = rbind(combDat, df)
    })
    
  }
  return(combDat)
}

ampseq_new = combDat(path = "ampseq_vivacity_found")
ampseq_new$Batch = "New"

#  make edited list of samples to exclude from summary of old samples
ampseq_new_samples = unique(ampseq_new$OrigName)
ampseq_new_samples = gsub("\\_S.*", "", ampseq_new_samples)
ampseq_new_samples = gsub("_", "-", ampseq_new_samples)

getOldSampleslist = function(targetPattern, targDir, newSamples) {
  combSamples = c()
  for ( curPattern in targetPattern) {
    dirsList = list.files(targDir, pattern = curPattern )
    for (curDir in dirsList) {
      curName = gsub("\\_S.*", "", curDir)
      curName = gsub("_", "-", curName)
      if (!curName%in%newSamples & !grepl("merge", curName, ignore.case = T) & !grepl("Water", curName, ignore.case = T))
        combSamples = c(combSamples, curDir)
    }
  }
  return(combSamples)
}

oldAmpList = getOldSampleslist(targetPattern=targetPattern, targDir="ampseq_vivacity_old", newSamples = ampseq_new_samples)

combOldDat = function(path, samplesList) {
  combDat = data.frame(matrix(nrow = 0, ncol =0))
  samples = list.files(path)
  samples = samples[samples%in%samplesList]
  for (curSample in samples) {
    try({
      curDf = paste0(path, "/", curSample, "/filtered.csv")
      df = read.csv(curDf)
      exactSample = gsub("_", "-", curSample)
      eaxactSample = sub(".*EHC", "EHC", curSample)
      df$OrigName = curSample
      df$ExactSample = exactSample
      sampNamelist = strsplit(exactSample, "-")
      sampleName =  sampNamelist[[1]][1:3]
      sampleName = paste(sampleName, collapse = "-")
      df$Sample = sampleName
      df$Samp_Pos_Ref_Alt = paste(df$Sample, df$POSITION, df$REF.NT, df$VAR.NT, sep = "__")
      
      combDat = rbind(combDat, df)
    })
    
  }
  return(combDat)
}

ampseq_old = combOldDat(path = "ampseq_vivacity_old", samplesList = oldAmpList)
ampseq_old$Batch = "Old"

# save ampseq
dir.create("snvs_comb_res")

write.csv(ampseq_new, "snvs_comb_res/ampseq_found.csv", row.names = F)
write.csv(ampseq_old, "snvs_comb_res/ampseq_old.csv", row.names = F)

## metaseq samples
metaseq_new = combDat(path = "metaseq_vivacity_found")
metaseq_new$Batch = "New"

#  make edited list of samples to exclude from summary of old samples
metaseq_new_samples = unique(metaseq_new$OrigName)
metaseq_new_samples = gsub("\\_S.*", "", metaseq_new_samples)
metaseq_new_samples = gsub("_", "-", metaseq_new_samples)

oldMetaList = getOldSampleslist(targetPattern=targetPattern, targDir="metaseq_vivacity_old", newSamples = metaseq_new_samples)

metaseq_old = combOldDat(path = "metaseq_vivacity_old", samplesList = oldMetaList)
metaseq_old$Batch = "Old"

write.csv(metaseq_new, "snvs_comb_res/metaseq_found.csv", row.names = F)
write.csv(metaseq_old, "snvs_comb_res/metaseq_old.csv", row.names = F)


