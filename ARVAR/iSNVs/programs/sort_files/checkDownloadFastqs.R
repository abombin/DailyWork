curDat = read.csv("SeqSamples_Eval_Depth10.csv")
# metaseq
missingSamples = unique(curDat$Missing_Metaseq)
missingSamples = missingSamples[!is.na(missingSamples)]

spreadList = function(x) {
  combList = character()
  for ( i in x ) {
    curList = strsplit(i, ";")[[1]]
    combList = c(combList, curList)
  } 
  return(combList)
}

curMisLibs = spreadList(missingSamples)

meta1  = list.files("metaseq_pass_fastqs")
meta2 = list.files("metaseq_test_samples", recursive = T)
newMeta = c(meta1, meta2)

missingPerLib = function(missingList, curFiles) {
  combMissing = character()
  for ( i in missingList ) {
    curSeqCh = strsplit(i, "_")[[1]]
    curMainStr = curSeqCh[1:3]
    curMainStr = paste(curMainStr, collapse = ".")
    if (length(curSeqCh) > 3) {
      curLib = curSeqCh[4]
      #curLib = gsub("[Ll]", "*",  curLib )
      curPattern = paste(curMainStr, curLib, sep = ".")
    } else {
      curPattern = paste0(curMainStr, "(?![_-]L)")
    }
    
    if (!any(grepl(curPattern, curFiles, ignore.case = T, perl = T))) {
      combMissing =  c(combMissing, i)
    }
  }
  return(combMissing)
}

grepl("EHC-C19_2236B(?![_-]L)", "EHC-C19_2236B_1_R1", ignore.case = T, perl = T)
missingLibs = missingPerLib(missingList = curMisLibs, curFiles = newMeta)
length(missingLibs)


missingPerSamp = function(missingList, curFiles) {
  combMissing = character()
  for ( i in missingList ) {
    curSeqCh = strsplit(i, "_")[[1]]
    curMainStr = curSeqCh[1:3]
    curMainStr = paste(curMainStr, collapse = ".")
    
    if (!any(grepl(curMainStr, curFiles, ignore.case = T, perl = T))) {
      combMissing =  c(combMissing, i)
    }
  }
  return(combMissing)
}

missSamples = missingPerSamp(missingList = curMisLibs, curFiles = newMeta)

# ampseq
ampSamples = unique(curDat$Missing_Ampseq)
missingSamples = ampSamples[!is.na(ampSamples)]
curMisLibs = spreadList(missingSamples)

#
amp1 = list.files("amp_fail_files")
amp2 = list.files("amp_miss_files")
amp3 = list.files("amp_pass_files")
amp4 = list.files("amp_spResolve_files")
newAmp = c(amp1, amp2, amp3, amp4)

missingLibs = missingPerLib(missingList = curMisLibs, curFiles = newAmp)
missSamples = missingPerSamp(missingList = curMisLibs, curFiles = newAmp)