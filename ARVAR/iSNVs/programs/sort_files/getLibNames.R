setwd("/home/ubuntu/extraVol/ARVAR/iSNVs")
ludy_df = read.csv("ludy_sequenced_samples.csv")

colnames(ludy_df)

getType = function(df, typeCol, depthCol, i, minDepth) {
  df[i, typeCol] = tolower(df[i, typeCol])
  df[i, depthCol] = as.numeric(as.character(df[i, depthCol]))
  if (df[i, typeCol] == "metaseq" && !is.na(df[i, depthCol]) && df[i, depthCol] > minDepth) {
    curVal = "metaseq"
  } else if (df[i, typeCol] == "ampseq" && !is.na(df[i, depthCol]) && df[i, depthCol] > minDepth) {
    curVal = "ampseq"
  } else {
    curVal = NA
  }
  return(curVal)
}

getLibNames = function(libList, curSample) {
  metaList = character()
  ampList = character()
  for ( j in 1:length(libList) ) {
    lib = libList[j]
    if ( names(libList)[j] == "L1") {
      newName = paste0(curSample)
    } else {
      newName = paste0(curSample, "_", names(libList)[j])
    }
    if (!is.na(lib) && lib == "metaseq") {
      metaList = c(metaList, newName)
    } else if (!is.na(lib) && lib == "ampseq") {
      ampList = c(ampList, newName)
    }
  }
  if (length(metaList) > 0) {
    metaList = metaList
  } else {
    metaList = NA
  }
  
  if (length(ampList) > 0) {
    ampList = ampList
  } else {
    ampList = NA
  }
  return(list(metaList, ampList))
}

sumTable = function(df, minDepth) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df) ) {
    
    curSample = df$Sample_ID[i]
    
    L1 = getType(df = df , typeCol = "sequencing.type", depthCol = "Depth", i = i, minDepth = minDepth)
    L2 = getType(df = df , typeCol = "sequencing.type2", depthCol = "Depth4", i = i, minDepth = minDepth)
    L3 = getType(df = df , typeCol = "Sequencing.type5", depthCol = "Depth7", i = i, minDepth = minDepth)
    L4 = getType(df = df , typeCol = "Sequencing.type8", depthCol = "Depth10", i = i, minDepth = minDepth)
    L5 = getType(df = df , typeCol = "Sequencing.type3", depthCol = "Depth2", i = i, minDepth = minDepth)
    LAmp = getType(df = df , typeCol = "Sequencing.type4", depthCol = "Depth3", i = i, minDepth = minDepth)
    
    libList = c(L1, L2, L3, L4, L5, LAmp)
    names(libList) = c('L1', 'L2', 'L3', 'L4', 'L5', 'LAmp')
    
   combLib = getLibNames(libList = libList, curSample = curSample)
   curMeta = paste(combLib[[1]], collapse = ";")
   curAmp = paste(combLib[[2]], collapse = ";")
   
   curDat = data.frame(curSample, curMeta, curAmp)
   combDat = rbind(combDat, curDat)
  }
  colnames(combDat) = c("Sample_id", "MetaseqNames", "AmpseqNames")
  return(combDat)
}

combDat = sumTable(df = ludy_df, minDepth = 10)

write.csv(combDat, "SeqSamplesNames_Depth10.csv", row.names = F)