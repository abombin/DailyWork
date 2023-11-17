library(dplyr)

fastqs = read.csv("metaseq_dx_fastqs.csv")

metaseq = read.csv("metaseq_libs_fail.csv")

exclList = c("EHC_C19_1407E_L2", "EHC_C19_1445Q_L2")
exclLibList = "EHC.C19.1407E.2|EHC.C19.1445Q.2"
metaseqFilt = metaseq[!metaseq$combSeqName%in%exclList,]
exclSamp = metaseq[metaseq$combSeqName%in%exclList,]

metaSpike = metaseqFilt[(metaseqFilt$Spike.in.confirmed == "Y"),]
metaSpFail = metaseqFilt[!(metaseqFilt$Spike.in.confirmed == "Y"),]
rownames(metaSpFail) = NULL

convertFileSize = function(fastqs) {
  fastqs$Size = NA
  for ( i in 1:nrow(fastqs) ) {
    if (fastqs$V5[i] == "GB") {
      fastqs$Size[i] = fastqs$V4[i] * 1073741824
    } else if (fastqs$V5[i] == "MB") {
      fastqs$Size[i] = fastqs$V4[i] * 1048576
    } else if (fastqs$V5[i] == "KB") {
      fastqs$Size[i] = fastqs$V4[i] * 1024
    } else if (fastqs$V5[i] == "bytes") {
      fastqs$Size[i] = fastqs$V4[i]
    }
  }
  return(fastqs)
}

fastqs = convertFileSize(fastqs = fastqs)

getBiggestFiles = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$V6%in%curFilesList,]
    size_mean <- curFastqs %>%
      group_by(V2) %>%
      summarise(SumValue = mean(Size))
    max_mean = max(size_mean$SumValue)
    max_df = size_mean[size_mean$SumValue ==  max_mean,]
    max_gr = max_df$V2[1]
    cur_fastqs_max = curFastqs[curFastqs$V2 ==  max_gr,]
    cur_fastqs_max$Lib_name = df$combSeqName[i]
    combDat = rbind(combDat, cur_fastqs_max)
  }
  return(combDat)
}

bigMeta = getBiggestFiles(df = metaSpike, fastqs = fastqs)
bmSum = data.frame(table(bigMeta$Lib_name))

#add excluded samples
exclFastq = fastqs[grepl(exclLibList, fastqs$V6),]
exclFastq$Lib_name = paste0(exclFastq$Sample_id, "_L2")
bigMetaComb = rbind(bigMeta, exclFastq)

#write.csv(bigMetaComb, "meta_lib_reolv_spike_size.csv", row.names = F)

# resolve samples that failed spike in
downloadFailed = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$V6%in%curFilesList,]
    curFastqs$Lib_name = df$combSeqName[i]
    curFastqs$NewDir = paste(curFastqs$V2, curFastqs$Lib_name, sep = "__")
    combDat = rbind(combDat, curFastqs)
    newDirs = unique(curFastqs$NewDir)
    for (newDir in newDirs) {
      curDir = paste0("metaseq_test_samples/", newDir, "/")
      dir.create(curDir, recursive = T, showWarnings = F)
      curFiles = curFastqs[curFastqs$NewDir == newDir,]$V6
      for (curFile in curFiles) {
        cmd_str = paste0("dx download ", curFile, " -o ", curDir)
        #system(cmd_str)
      }
    }
  }
  return(combDat)
}

failedFastqs = downloadFailed(df = metaSpFail, fastqs = fastqs)

# can be parallel
statFailedLib = function(inDir) {
  dirsList = list.files(inDir)
  outDir = "meta_fail_stats/"
  dir.create(outDir, recursive = T, showWarnings = F)
  for ( curDir in dirsList) {
    cmd_str = paste0("seqkit stats -a ", inDir, "/", curDir, "/*")
    curStats = system(cmd_str, intern = T)
    outFile = paste0(outDir, curDir)
    write.table(curStats, outFile, row.names = F, col.names = F, quote = F)
  }
}

statFailedLib(inDir = "metaseq_test_samples")

combStats = function(inDir) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  filesList = list.files(inDir)
  for ( curFile in filesList ) {
    inFile = paste0(inDir, "/", curFile)
    curDf = read.table(inFile, T)
    curDf$NewDir = curFile
    combDat = rbind(combDat, curDf)
  }
  return(combDat)
}

statsMetaFailed = combStats(inDir = "meta_fail_stats")

failedFastqs$FileID = basename(failedFastqs$V6)
failedFastqs$FullID = paste(failedFastqs$NewDir, failedFastqs$FileID, sep = "__")
statsMetaFailed$FileID = basename(statsMetaFailed$file)
statsMetaFailed$FullID = paste(statsMetaFailed$NewDir, statsMetaFailed$FileID, sep = "__")

failedFastqsStats = unique(plyr::join(failedFastqs, statsMetaFailed, by = "FullID", type = "left", match = 'all'))

failedFastqsStats$IQR = failedFastqsStats$Q3 - failedFastqsStats$Q1

statsFilt = failedFastqsStats[, c("V2",  "V3", "V4", "V5", "V6", "Size", "Lib_name", "min_len", "avg_len", "max_len", "Q1",  "Q2", "Q3" ,"IQR" )]

colnames(statsFilt)[1:5] = c("Dx_date", "Dx_time", "Dx_size", "Dx_unit", "Dx_path")

write.csv(statsFilt, "metaseq_failed_stats.csv", row.names = F)

system("aws s3 cp metaseq_failed_stats.csv s3://abombin/ARVAR/iSNVs/")