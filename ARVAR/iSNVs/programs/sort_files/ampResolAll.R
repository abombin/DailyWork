library(dplyr)

fastqs = read.csv("ampseq_combined_paths.csv")
fastqs$Pbnas_path = gsub('\\\\', "/", fastqs$Pbnas_path)
fastqs = fastqs[grepl("fastq.gz", fastqs$Pbnas_path),]

ampseq = read.csv("ampseq_libs_fail.csv")


ampSpike = ampseq[(ampseq$Spike.in.confirmed == "AMPSEQ"),]
ampSpFail = ampseq[!(ampseq$Spike.in.confirmed == "AMPSEQ"),]
rownames(ampSpFail) = NULL

getBiggestFiles = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$Pbnas_path%in%curFilesList,]
    curFastqs$Group = NA
    for (j in 1:nrow(curFastqs)) {
      pathList = strsplit(curFastqs$Pbnas_path[j], "/")[[1]]
      grInd = length(pathList) - 1
      elem =  pathList[grInd]
      curFastqs$Group[j] =  elem
    }

    
    size_mean <- curFastqs %>%
      group_by(Group) %>%
      summarise(SumValue = mean(File_size))
    max_mean = max(size_mean$SumValue)
    max_df = size_mean[size_mean$SumValue ==  max_mean,]
    max_gr = max_df$Group[1]
    cur_fastqs_max = curFastqs[curFastqs$Group ==  max_gr,]
    cur_fastqs_max$Lib_name = df$combSeqName[i]
    combDat = rbind(combDat, cur_fastqs_max)
  }
  return(combDat)
}

bigamp = getBiggestFiles(df = ampSpike, fastqs = fastqs)
bmSum = data.frame(table(bigamp$Lib_name))
write.csv(bigamp, "amp_fail_spResolv.csv", row.names = F)
system("aws s3 cp amp_fail_spResolv.csv s3://abombin/ARVAR/iSNVs/")


bigampFail = unique(getBiggestFiles(df = ampSpFail, fastqs = fastqs))


splitFailFiles = function(df, fastqs) {
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

failedFilesFastqs = unique(splitFailFiles(df = ampSpFail, fastqs = fastqs))
rownames(failedFilesFastqs) = NULL

# write.csv(failedFilesFastqs, "ampseq_failed_diffSize.csv", row.names = F)
# system("aws s3 cp ampseq_failed_diffSize.csv s3://abombin/ARVAR/iSNVs/")

# resolve missing
ampMiss = read.csv("ampseq_libs_miss.csv")

splitMissFiles = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$Pbnas_path%in%curFilesList,]
    curFastqs$Group = NA
    curFastqs$combSeqName = NA
    curFastqs$combOrigMiss = NA
    curFastqs$Spike.in.confirmed = NA
    
    curFastqs$combSeqName = df$combSeqName[i]
    curFastqs$combOrigMiss = df$combOrigMiss[i]
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


missFastqs = unique(splitMissFiles(df = ampMiss , fastqs = fastqs))

# get unique miss files
getUniqueFastqs = function(df) {
  uniqueList = character()
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df) ) {
    curRow = df[i ,]
    curFile = paste(basename(df$Pbnas_path[i]), df$File_size[i], sep = "__")
    if (!curFile %in% uniqueList) {
      combDat = rbind(combDat, curRow )
      uniqueList = c(uniqueList, curFile)
    }
  }
  return(combDat)
}

uniqueMiss = getUniqueFastqs(df = missFastqs)

sum(uniqueMiss$File_size) / 1000

# get unique failed samples
uniqueFail = getUniqueFastqs(df = failedFilesFastqs)

# write.csv(uniqueMiss, "amp_uniqueMiss.csv", row.names = F)
# write.csv(uniqueFail, "amp_uniqueFail.csv", row.names = F)
# 
# system("aws s3 cp amp_uniqueMiss.csv s3://abombin/ARVAR/iSNVs/")
# system("aws s3 cp amp_uniqueFail.csv s3://abombin/ARVAR/iSNVs/")

filterUniqueFail = function(df) {
  libNames = unique(df$combSeqName)
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( libName in libNames ) {
    curDf = df[df$combSeqName == libName,]
    if (nrow(curDf) > 2) {
      curSelDf = curDf[grepl("may222023", curDf$Group),]
      if (nrow(curSelDf) < 2 ) {
        curSelDf = curDf
      }
    } else {
      curSelDf = curDf
    }
    combDat = rbind(combDat, curSelDf)
  }
  return(combDat)
}

uniqueFailFiltr = filterUniqueFail(df = uniqueFail)

write.csv(uniqueFailFiltr, "amp_uniqueFailFiltr.csv", row.names = F)
system("aws s3 cp amp_uniqueFailFiltr.csv s3://abombin/ARVAR/iSNVs/")
