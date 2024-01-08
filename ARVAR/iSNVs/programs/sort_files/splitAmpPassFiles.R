ampPass = read.csv("~/extraVol/ARVAR/iSNVs/ampseq_libs_pass.csv")

fastqs = read.csv("ampseq_combined_paths.csv")
fastqs$Pbnas_path = gsub('\\\\', "/", fastqs$Pbnas_path)
fastqs = fastqs[grepl("fastq.gz", fastqs$Pbnas_path),]

splitPassFiles = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$Pbnas_path%in%curFilesList,]
    curFastqs$Group = NA
    curFastqs$combSeqName = NA
    curFastqs$Spike.in.confirmed = NA
    
    curFastqs$combSeqName = df$combSeqName[i]
    curFastqs$Spike.in.confirmed = df$Spike.in.confirmed[i]
    for (j in 1:nrow(curFastqs)) {
      pathList = strsplit(curFastqs$Pbnas_path[j], "/")[[1]]
      grInd = length(pathList) - 1
      elem =  pathList[grInd]
      curFastqs$Group[j] =  elem
    }
    combDat = rbind(combDat, curFastqs)
  }
  return(combDat)
}

passFastqs = unique(splitPassFiles(df = ampPass, fastqs = fastqs))
write.csv(passFastqs, "ampseq_pass_files.csv", row.names = F)

system("aws s3 cp ampseq_pass_files.csv s3://abombin/ARVAR/iSNVs/")