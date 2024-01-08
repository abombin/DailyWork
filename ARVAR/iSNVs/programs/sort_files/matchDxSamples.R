

combAll = function(filesList, inDir) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for (curFile in filesList) {
    inFile = paste0(inDir, "/", curFile)
    curDf = read.table(inFile)
    curDf$Sample_id = curFile
    combDat = rbind(combDat, curDf)
  }
  return(combDat)
}

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

# add dx fastqs data
addDx = function(mainDf, dxDf) {
  mainDf$Meta_dx = NA
  mainDf$Meta_dx_fastq = NA

  for ( i in 1:nrow(mainDf) ) {
    curSample = mainDf$Sample_id[i]
    if (!is.na(mainDf$Missing_Metaseq[i])) {
      curDxDf = dxDf[dxDf$Sample_id == curSample,]
      if (nrow(curDxDf) > 0) {
        curFiles = curDxDf$V6
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
  }
  return(mainDf)
}

# dxMissing = function(df) {
#   df$Missing_meta_dx = NA
#   for ( i in 1:nrow(df) ) {
#     metaseqStr = df$Missing_Metaseq[i]
#     if (!is.na(metaseqStr)) {
#       curMetaseq = strsplit(metaseqStr, ";")[[1]]
#       curDx = strsplit(df$Meta_dx[i], ";")[[1]]
#       curMissingMetaseq = character()
#       for (curSeq in curMetaseq) {
#         curSeqCh = strsplit(curSeq, "_")[[1]]
#         curMainStr = curSeqCh[1:3]
#         curMainStr = paste(curMainStr, collapse = ".")
#         if (length(curSeqCh) > 3) {
#           curLib = curSeqCh[4]
#           curLib = gsub("[Ll]", "*",  curLib )
#           curPattern = paste(curMainStr, curLib, sep = ".")
#         } else {
#           curPattern = curMainStr
#         }
# 
#         if (!any(grepl(curPattern, curDx))) {
#           curMissingMetaseq =  c(curMissingMetaseq, curSeq)
#         }
#         
#         if (length(curMissingMetaseq) == 0) {
#           curMissingMetaseq = NA
#         }
#         if (!is.na(curMissingMetaseq)) {
#           curMissingMetaseq = paste(curMissingMetaseq, collapse = ";")
#         }
#         df$Missing_meta_dx[i] = curMissingMetaseq
#       }
#     }
#   }
#   return(df)
# }

dxMissing = function(df) {
  df$Missing_meta_dx = NA
  for ( i in 1:nrow(df) ) {
    metaseqStr = df$Missing_Metaseq[i]
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
    metaseqStr = df$Missing_Metaseq[i]
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
    metaseqStr = df$Missing_Metaseq[i]
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
    if ( length(curFastq) == 2 & any(grepl("_R1", curFastq)) & any(grepl("_R2", curFastq)) & df$Spike.in.confirmed[i] == "Y" ) {
      combDat = rbind(combDat, df[i,])
    } else {
      failDat = rbind(failDat, df[i,])
    }
  }
  dataList = list(passDat = combDat, failDat = failDat)
  return(dataList)
}

adjustL2 = function(df) {
  df$Original_miss = df$Missing_Metaseq
  for ( i in 1:nrow(df) ) {
    if (!grepl("_L2", df$Missing_Metaseq[i])) {
      df$Missing_Metaseq[i] = paste0(df$Missing_Metaseq[i], "_L2")
    }
  }
  return(df)
}

# current samples table
curDat = read.csv("SeqSamples_Eval_Depth10.csv")
metaSamples = curDat[!is.na(curDat$MetaseqNames),]
metaSamples = metaSamples[, 1:4]
rownames(metaSamples) = 1:nrow(metaSamples)

# get metaseq data from dx
filesList = list.files("dx_metaseq_paths")
combDat = combAll(filesList = filesList , inDir = "dx_metaseq_paths")
fastqs = combDat[grepl("fastq.gz", combDat$V6),]

#write.csv(fastqs, "metaseq_dx_fastqs.csv", row.names = F)

# update main table with data from dx
metaDxUpdate = addDx(mainDf = metaSamples, dxDf = fastqs)
rownames(metaDxUpdate) = 1:nrow(metaDxUpdate)

# update the table if anything is missing after dx search
metaDxUpdate = dxMissing(df = metaDxUpdate)

# check if any needed samples are not in DNA nexus
missingDx = metaDxUpdate[!is.na(metaDxUpdate$Missing_meta_dx),]

# extract file names that need to be downloaded
metaseqDxSamples = metaDxUpdate[!is.na(metaDxUpdate$Missing_Metaseq),]
metaseqDxSamples = metaseqDxSamples[is.na(metaseqDxSamples$Missing_meta_dx),]
rownames(metaseqDxSamples) = NULL
metaseqDxDownload = selectPathsMetaseq(df = metaseqDxSamples)

# represent as one row per Lib with download samples
metaseqDownloadLibs = getFastqPerLib(df = metaseqDxDownload )
write.csv(metaseqDownloadLibs, "metaseq_additional_samples.csv", row.names = F)
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

missingDx$Meta_dx_download = missingDx$Meta_dx_fastq
missingDx = adjustL2(df = missingDx)
missDownloadLibs = getFastqPerLib(df = missingDx)
missDownloadLibs$combSeqName[3:nrow(missDownloadLibs)] = gsub("_L2", "", missDownloadLibs$combSeqName[3:nrow(missDownloadLibs)])

missDLSpike = plyr::join(missDownloadLibs, postseq_sel, by = "combSeqName", type = "left", match = "all")

metaMissCheck = getClearDownloadSamples(df = missDLSpike)
metaMissPass = metaMissCheck[[1]]
metaMissFail = metaMissCheck[[2]]
missingDx = metaDxUpdate[!is.na(metaDxUpdate$Missing_meta_dx),]
metaMissFail$combFastq[1] = missingDx$Meta_dx_fastq[1]
metaMissFail$combFastq[2] = missingDx$Meta_dx_fastq[2]

# combine automatically resolved and those that have to be reviewd further
combPassMetaseq = rbind(metaPass, metaMissPass)
combFailMeta = rbind(metaFail, metaMissFail)

nrow(combPassMetaseq)
nrow(combFailMeta)

write.csv(combPassMetaseq, "metaseq_libs_pass.csv", row.names = F)
write.csv(combFailMeta, "metaseq_libs_fail.csv", row.names = F)

system("aws s3 cp metaseq_libs_fail.csv s3://abombin/ARVAR/iSNVs/")
system("aws s3 cp metaseq_libs_pass.csv s3://abombin/ARVAR/iSNVs/")