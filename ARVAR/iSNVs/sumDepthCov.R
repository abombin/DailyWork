library(doParallel)
library(foreach)


getFilesList = function(inPath, pattern) {
  filesList = list.files(inPath, pattern = pattern, full.names = T, recursive = T)
  return(filesList)
}


getStats = function(inFile) {
  cmd_str = paste0("/home/ubuntu/extraVol/bio-soft/samtools coverage --depth 0 ", inFile)
  statsStr = system(cmd_str, intern = T)
  outStrs = strsplit(statsStr, " ")
  curStats = strsplit(outStrs[[2]], "\t")[[1]]
  curCov = curStats[6]
  curDepth = curStats[7]
  #curFile = basename(inFile)
  curFile = gsub("\\/output.bam", "", inFile)
  curFile = basename(curFile)
  df = data.frame(curFile, curCov, curDepth)
  colnames(df) = c("Sample_id", "Coverage", "Mean_depth")
  return(df)
}


getAllStats = function(inFile) {
  cmd_str = paste0("/home/ubuntu/extraVol/bio-soft/samtools coverage --depth 0 ", inFile)
  statsStr = system(cmd_str, intern = T)
  outStrs = strsplit(statsStr, " ")
  curStats = strsplit(outStrs[[2]], "\t")[[1]]
  curNames = strsplit(outStrs[[1]], "\t")[[1]]
  #curFile = basename(inFile)

  df = data.frame(t(curStats))
  colnames(df) = curNames
  colnames(df)[1] = "OrigSamp"
  return(df)
}

# run parallel
runPar = function(filesList) {
  useCores = 30
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)
  
  combDat = foreach(i = filesList, .combine  = 'rbind') %dopar% {
    getAllStats(inFile = i)
  }
  parallel::stopCluster(cl = cl)
  return(combDat)
}

#Ampseq overlap
# ampFiles = getFilesList(inPath="IntraSnv_ampseq_overlap", pattern="output.bam$")
# ampStats = runPar(filesList = ampFiles)
# write.csv(ampStats, "IntraSnv_ampseq_overlap/ampseq_stats_2.csv", row.names = F)

#metaseq overlap
# metaFiles = getFilesList(inPath="IntraSnv_metaseq_overlap", pattern="output.bam$")
#metaSel = metaFiles[grepl("EHC-C19-1636Z", metaFiles)]
# metaStats = runPar(filesList = metaFiles)
# write.csv(metaStats, "IntraSnv_metaseq_overlap/metaseq_stats_2.csv", row.names = F)



#Ampseq overlap
ampFiles = getFilesList(inPath="IntraSnv_ampseq_all", pattern="output.bam$")
ampStats = runPar(filesList = ampFiles)
write.csv(ampStats, "IntraSnv_results/ampseq_stats_2.csv", row.names = F)

#metaseq overlap
metaFiles = getFilesList(inPath="IntraSnv_metaseq_all", pattern="output.bam$")
metaStats = runPar(filesList = metaFiles)
write.csv(metaStats, "IntraSnv_results/metaseq_stats_2.csv", row.names = F)
