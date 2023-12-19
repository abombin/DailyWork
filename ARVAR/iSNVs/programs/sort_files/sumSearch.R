curDat = read.csv("SeqSamples_Eval_Depth10.csv")
selDat = curDat[, 1:4]

meta1  = list.files("metaseq_pass_fastqs")
meta2 = list.files("metaseq_test_samples", recursive = T)
newMeta = c(meta1, meta2)

editName = function(filesList) {
  combList = character()
  for ( i in filesList ) {
    newName = gsub("^.*\\/", "", i)
    newName = gsub("-", "_", newName)
    newName = gsub('\\_S.*', "", newName)
    combList = c(combList, newName)
  }
  unList = unique(combList)
  return(unList)
}


editMetaList = function(df, filesList) {
  df$Found_Samples = NA
  for ( i in 1:nrow(df) ) {
    presSamples = strsplit(df$MetaFound[i], ";")[[1]]
    presPattern = gsub("_", ".", df$Sample_id[i])
    foundSamples = filesList[grepl(presPattern, filesList)]
    newSamples = foundSamples[!(foundSamples%in%presSamples)]
    if (length(newSamples) > 0) {
      newSampStr = paste(newSamples, collapse = ";")
    } else {
      newSampStr = NA
    }
    df$Found_Samples[i] = newSampStr 
  }
  return(df)
}

curMeta = editName(filesList = newMeta)
editMetaDf = editMetaList(df = selDat, filesList = curMeta)
colnames(editMetaDf)

metaAllSamp = editMetaDf[, c("Sample_id", "MetaseqNames", "MetaFound" , "Found_Samples")]
colnames(metaAllSamp) = c("Sample_id", "Ludy_Lib_Names", "Orig_Lib_Found", "New_Lib_Found")

# write.csv(metaAllSamp, "metaseq_search_result.csv", row.names = F)
# system("aws s3 cp metaseq_search_result.csv s3://abombin/ARVAR/iSNVs/")


## AMPSEQ

# need to upload new files and make a new check

curDat = read.csv("SeqSamples_Eval_Depth10.csv")
selDat = curDat[, c(1,5:7)]

amp1 = list.files("amp_fail_files")
amp2 = list.files("amp_miss_files")
amp3 = list.files("amp_pass_files")
amp4 = list.files("amp_spResolve_files")
amp5 = list.files("ampseq_additional")
newAmp = c(amp1, amp2, amp3, amp4, amp5)

editName = function(filesList) {
  combList = character()
  for ( i in filesList ) {
    newName = gsub("^.*\\/", "", i)
    newName = gsub("-", "_", newName)
    newName = gsub('\\_S.*', "", newName)
    combList = c(combList, newName)
  }
  unList = unique(combList)
  return(unList)
}


editAmpList = function(df, filesList) {
  df$Found_Samples = NA
  for ( i in 1:nrow(df) ) {
    presSamples = strsplit(df$AmpFound[i], ";")[[1]]
    presPattern = gsub("_", ".", df$Sample_id[i])
    foundSamples = filesList[grepl(presPattern, filesList)]
    newSamples = foundSamples[!(foundSamples%in%presSamples)]
    if (length(newSamples) > 0) {
      newSampStr = paste(newSamples, collapse = ";")
    } else {
      newSampStr = NA
    }
    df$Found_Samples[i] = newSampStr 
  }
  return(df)
}

curAmp = editName(filesList = newAmp)
editAmpList = editAmpList(df = selDat, filesList = curAmp)

ampAllSamp = editAmpList[, c("Sample_id", "AmpseqNames", "AmpFound" , "Found_Samples")]
colnames(ampAllSamp) = c("Sample_id", "Ludy_Lib_Names", "Orig_Lib_Found", "New_Lib_Found")

 write.csv(ampAllSamp, "ampseq_search_result.csv", row.names = F)
 system("aws s3 cp ampseq_search_result.csv s3://abombin/ARVAR/iSNVs/")