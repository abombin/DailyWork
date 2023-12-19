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

curMeta = editName(filesList = newMeta)

additionalDat = read.csv("metaseq_additional_samples.csv")

toDownload = additionalDat[!additionalDat$combSeqName%in%curMeta,]

for ( i in 1:nrow(toDownload) ) {
  curFiles = strsplit(toDownload$combFilesStr[i], ";")[[1]]
  for (curFile in curFiles) {
    if ( !grepl("Ludy_Apr242023", curFile, ignore.case = T) & !grepl("for_andre", curFile, ignore.case = T) ) {
      cmd_str = paste0("dx download ", curFile, " -f -o  metaseq_pass_fastqs/")
      try({
        system(cmd_str)
      })
    }

  }
}

## AMPSEQ
amp1 = list.files("amp_fail_files")
amp2 = list.files("amp_miss_files")
amp3 = list.files("amp_pass_files")
amp4 = list.files("amp_spResolve_files")
newAmp = c(amp1, amp2, amp3, amp4)

curAmp = editName(filesList = newAmp)

additionalDat = read.csv("ampseq_additional_samples.csv")

toDownload = additionalDat[!additionalDat$combSeqName%in%curAmp,]
rownames(toDownload) = NULL