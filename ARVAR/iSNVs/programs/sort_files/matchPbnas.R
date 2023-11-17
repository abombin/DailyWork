fastqs = read.csv("ampseq_combined_paths.csv")
fastqs$Pbnas_path = gsub('\\\\', "/", fastqs$Pbnas_path)
fastqs = fastqs[grepl("fastq.gz", fastqs$Pbnas_path),]

# current samples table
curDat = read.csv("SeqSamples_Eval_Depth10.csv")
ampSamples = curDat[!is.na(curDat$Missing_Ampseq),]
ampSamples = ampSamples[, c(1, 5:7)]
rownames(ampSamples) = 1:nrow(ampSamples)


formatNames = function(samplesList) {
  newNames = character()
  samplesList = basename(samplesList)
  sampleNames = gsub("-", "_", samplesList)
  sampleNames = sub(".*EHC", "EHC", sampleNames)
  for ( curSample in sampleNames ) {
    sampNamelist = strsplit(curSample, "_")
    sampleName =  sampNamelist[[1]][1:4]
    sampleName = paste(sampleName, collapse = "_")
    sampleName = gsub('\\_S.*', "", sampleName)
    if (grepl("WATER", sampleName, ignore.case = T) == F) {
      newNames = c(newNames, sampleName)
    }
  }
  return(newNames)
}

addDx = function(mainDf, dxDf) {
  mainDf$Meta_dx = NA
  mainDf$Meta_dx_fastq = NA
  
  for ( i in 1:nrow(mainDf) ) {
    curSample = mainDf$Sample_id[i]
      curDxDf = dxDf[dxDf$Sample_id == curSample,]
      if (nrow(curDxDf) > 0) {
        curFiles = curDxDf$Pbnas_path
        curFilesString = paste(curFiles, collapse = ";")
        curNames = unique(formatNames(samplesList = curFiles))
        curNamesStr = paste(curNames, collapse = ";")
        
      } else {
        curStrings = NA
        curNamesStr = NA
      }
      mainDf$Meta_dx[i] = curNamesStr
      mainDf$Meta_dx_fastq[i] = curFilesString
  }
  return(mainDf)
}

dxMissing = function(df) {
  df$Missing_meta_dx = NA
  for ( i in 1:nrow(df) ) {
    metaseqStr = df$Missing_Ampseq[i]
    if (!is.na(metaseqStr)) {
      curMetaseq = strsplit(metaseqStr, ";")[[1]]
      curDx = strsplit(df$Meta_dx[i], ";")[[1]]
      curMissingMetaseq = character()
      for (curSeq in curMetaseq) {
        curSeqCh = strsplit(curSeq, "_")[[1]]
        curMainStr = curSeqCh[1:3]
        curMainStr = paste(curMainStr, collapse = ".")
        if (length(curSeqCh) > 3) {
          curLib = curSeqCh[4]
          #curLib = gsub("[Ll]", "*",  curLib )
          curPattern = paste(curMainStr, curLib, sep = ".")
        } else {
          curPattern = paste0(curMainStr, "(?![_-].*)")
        }
        
        if (!any(grepl(curPattern, curDx, ignore.case = T, perl = T))) {
          curMissingMetaseq =  c(curMissingMetaseq, curSeq)
        }
        
      }
      
      if (length(curMissingMetaseq) == 0) {
        curMissingMetaseq = NA
      }
      if (!is.na(curMissingMetaseq)) {
        curMissingMetaseq = paste(curMissingMetaseq, collapse = ";")
      }
      df$Missing_meta_dx[i] = curMissingMetaseq
    }
  }
  return(df)
}



selectPathsMetaseq = function(df) {
  df$Meta_dx_download = NA
  for ( i in 1:nrow(df)) {
    metaseqStr = df$Missing_Ampseq[i]
    if (!is.na(metaseqStr)) {
      curMetaseq = strsplit(metaseqStr, ";")[[1]]
      curDx = strsplit(df$Meta_dx_fastq[i], ";")[[1]]
      filesComb = character()
      for (curSeq in curMetaseq) {
        curSeqCh = strsplit(curSeq, "_")[[1]]
        curMainStr = curSeqCh[1:3]
        curMainStr = paste(curMainStr, collapse = ".")
        if (length(curSeqCh) > 3) {
          curLib = curSeqCh[4]
          #curLib = gsub("[Ll]", "*",  curLib )
          curLib = paste0(curLib, "[_-]S")
          curPattern = paste(curMainStr, curLib, sep = ".")
          curFiles = unique(curDx[grepl(curPattern, curDx, ignore.case = T)])
        } else if (length(curSeqCh) == 3) {
          curPattern = paste0(curMainStr, "[_-].*[_-]S")
          curFiles = unique(curDx[!grepl(curPattern, curDx, ignore.case = T)])
        }
        filesComb = c(filesComb, curFiles)
      }
      curFilesStr = paste(filesComb, collapse = ";")
      df$Meta_dx_download[i] = curFilesStr
    }
  }
  return(df)
}

getFastqPerLib = function(df) {
  combSeqName = character()
  combFastq = character()
  combFilesStr = character()
  for ( i in 1:nrow(df)) {
    metaseqStr = df$Missing_Ampseq[i]
    curMetaseq = strsplit(metaseqStr, ";")[[1]]
    curDx = strsplit(df$Meta_dx_download[i], ";")[[1]]
    curFastq = character()
    for (curSeq in curMetaseq) {
      curSeqCh = strsplit(curSeq, "_")[[1]]
      curMainStr = curSeqCh[1:3]
      curMainStr = paste(curMainStr, collapse = ".")
      if (length(curSeqCh) > 3) {
        curLib = curSeqCh[4]
        #curLib = gsub("[Ll]", "*",  curLib )
        curPattern = paste(curMainStr, curLib, sep = ".")
        curFiles = unique(curDx[grepl(curPattern, curDx, ignore.case = T)])
      } else if (length(curSeqCh) == 3) {
        curPattern = paste0(curMainStr, "[_-].*[_-]S")
        curFiles = unique(curDx[!grepl(curPattern, curDx, ignore.case = T)])
      }
      combSeqName = c(combSeqName, curSeq)
      FilesStr = paste( curFiles, collapse = ";")
      combFilesStr = c(combFilesStr, FilesStr)
    }
  }
  combData = data.frame(combSeqName, combFilesStr)
  return(combData)
}

getClearDownloadSamples = function(df) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  failDat = data.frame(matrix(nrow = 0, ncol = 0))
  df$Spike.in.confirmed = toupper(df$Spike.in.confirmed)
  df$Spike.in.confirmed[is.na(df$Spike.in.confirmed)] = 0
  df$Spike.in.confirmed = gsub(" ", "", df$Spike.in.confirmed)
  for ( i in 1:nrow(df) ) {
    curFastq = strsplit(df$combFastq[i], ";")[[1]]
    if ( length(curFastq) == 2 & any(grepl("_R1", curFastq)) & any(grepl("_R2", curFastq)) & df$Spike.in.confirmed[i] == "AMPSEQ" ) {
      combDat = rbind(combDat, df[i,])
    } else {
      failDat = rbind(failDat, df[i,])
    }
  }
  dataList = list(passDat = combDat, failDat = failDat)
  return(dataList)
}

adjustL2 = function(df) {
  df$Original_miss = df$Missing_Ampseq
  for ( i in 1:nrow(df) ) {
    if (!grepl("_L2", df$Missing_Ampseq[i])) {
      df$Missing_Ampseq[i] = paste0(df$Missing_Ampseq[i], "_L2")
    }
  }
  return(df)
}

ampUpdate = addDx(mainDf = ampSamples, dxDf = fastqs)

# update the table if anything is missing after dx search
metaDxUpdate = dxMissing(df = ampUpdate)

# check if any needed samples are not in DNA nexus
missingDx = metaDxUpdate[!is.na(metaDxUpdate$Missing_meta_dx),]

# extract file names that need to be downloaded
metaseqDxSamples = metaDxUpdate[!is.na(metaDxUpdate$Missing_Ampseq),]
metaseqDxSamples = metaseqDxSamples[is.na(metaseqDxSamples$Missing_meta_dx),]
rownames(metaseqDxSamples) = NULL
metaseqDxDownload = selectPathsMetaseq(df = metaseqDxSamples)

# represent as one row per Lib with download samples
metaseqDownloadLibs = getFastqPerLib(df = metaseqDxDownload )
write.csv(metaseqDownloadLibs, "ampseq_additional_samples.csv", row.names = F)

# check against Annes sheet
postseq = read.csv("Anne_postseq.csv")
postseq$Sample.ID = gsub("-", "_", postseq$Sample.ID)
postseq_sel = unique(postseq[, c("Sample.ID" , "Spike.in.confirmed" )])
colnames(postseq_sel)[1] = "combSeqName"
metaseqDLSpike = plyr::join(metaseqDownloadLibs, postseq_sel, by = "combSeqName", type = "left", match = "all")

metaCheck = getClearDownloadSamples(df = metaseqDLSpike)
metaPass = metaCheck[[1]]
metaFail = metaCheck[[2]]

# adjust for missingDX
missingDx = metaDxUpdate[!is.na(metaDxUpdate$Missing_meta_dx),]
rownames(missingDx) = NULL

# get new samples from missing
getNewSample = function(df) {
  df$newPbnas = NA
  for ( i in 1:nrow(df) ) {
    curSamples =  strsplit(df$Meta_dx[i], ";")[[1]]
    foundSamples = strsplit(df$AmpFound[i], ";")[[1]]
    newSamples = curSamples[!curSamples%in%foundSamples]
    if (length(newSamples) > 0) {
      newStr = paste(newSamples, collapse = ";")
      df$newPbnas[i] =  newStr
    } else {
      df$newPbnas[i] =  NA
    }

  }
  return(df)
}

editMissing = getNewSample(df = missingDx)

missingNew = editMissing[!is.na(editMissing$newPbnas),]
rownames(missingNew) = NULL
# save totally lost samples
ampTotLost =  editMissing[is.na(editMissing$newPbnas),]


# split libraries for missing
getFastqPerLibMiss = function(df) {
  combSeqName = character()
  combFastq = character()
  combOrigMiss = character()
  for ( i in 1:nrow(df)) {
    metaseqStr = df$newPbnas[i]
    curMetaseq = strsplit(metaseqStr, ";")[[1]]
    curDx = strsplit(df$Meta_dx_fastq[i], ";")[[1]]
    curOrigMiss = df$Missing_Ampseq[i]
    for (curSeq in curMetaseq) {
      curFastq = character()
      curSeqCh = strsplit(curSeq, "_")[[1]]
      curMainStr = curSeqCh[1:3]
      curMainStr = paste(curMainStr, collapse = ".")
      if (length(curSeqCh) > 3) {
        curLib = curSeqCh[4]
        #curLib = gsub("[Ll]", "*",  curLib )
        curPattern = paste(curMainStr, curLib, sep = ".")
        curFiles = unique(curDx[grepl(curPattern, curDx, ignore.case = T)])
      } else if (length(curSeqCh) == 3) {
        curPattern = paste0(curMainStr, "[_-].*[_-]S")
        curFiles = unique(curDx[!grepl(curPattern, curDx, ignore.case = T)])
      }
      curFastq = c(curFastq, curFiles)
      combFilesStr = paste(curFastq, collapse = ";")
    }
    combSeqName = c(combSeqName, curSeq)
    combFastq = c(combFastq, combFilesStr)
    combOrigMiss = c(combOrigMiss, curOrigMiss)
  }
  combData = data.frame(combSeqName, combFastq, combOrigMiss)
  return(combData)
}

misNewSplit  = getFastqPerLibMiss(df = missingNew)
misNewSpike = plyr::join(misNewSplit, postseq_sel, by = "combSeqName", type = "left", match = "all")
misNewSpike$Spike.in.confirmed = toupper(misNewSpike$Spike.in.confirmed)
#
newSpikeAmp = misNewSpike[misNewSpike$Spike.in.confirmed == "AMPSEQ",]
newSpikeAmp = newSpikeAmp[, c("combSeqName", "combFastq", "Spike.in.confirmed")]

finalAmpFail = rbind(metaFail, newSpikeAmp)
finalMiss = misNewSpike[!misNewSpike$Spike.in.confirmed == "AMPSEQ",]

write.csv(metaPass, "ampseq_libs_pass.csv", row.names = F)
write.csv(finalAmpFail, "ampseq_libs_fail.csv", row.names = F)
write.csv(finalMiss, "ampseq_libs_miss.csv", row.names = F)
write.csv(ampTotLost, "ampseq_libs_TotLost.csv", row.names = F)

