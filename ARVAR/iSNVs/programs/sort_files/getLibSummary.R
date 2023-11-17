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
    curVal = NULL
  }
  return(curVal)
}

sumTable = function(df, minDepth) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df) ) {
    
    curSample = df$Sample_ID[i]
    
    l1 = getType(df = df , typeCol = "sequencing.type", depthCol = "Depth", i = i, minDepth = minDepth)
    l2 = getType(df = df , typeCol = "sequencing.type2", depthCol = "Depth4", i = i, minDepth = minDepth)
    l3 = getType(df = df , typeCol = "Sequencing.type5", depthCol = "Depth7", i = i, minDepth = minDepth)
    l4 = getType(df = df , typeCol = "Sequencing.type8", depthCol = "Depth10", i = i, minDepth = minDepth)
    l5 = getType(df = df , typeCol = "Sequencing.type3", depthCol = "Depth2", i = i, minDepth = minDepth)
    l6 = getType(df = df , typeCol = "Sequencing.type4", depthCol = "Depth3", i = i, minDepth = minDepth)
   combLib = c(l1,l2,l3,l4,l5,l6)
   combLibSum = table(combLib)
   if ("metaseq"%in%names(combLibSum)) {
     curMetaseqCount = as.numeric(combLibSum["metaseq"])
   } else {
     curMetaseqCount = 0
   }
   if ("ampseq"%in%names(combLibSum)) {
     curAmpseqCount = as.numeric(combLibSum["ampseq"])
   } else {
     curAmpseqCount = 0
   }
   
   curData = data.frame(curSample, curMetaseqCount, curAmpseqCount)
   combDat = rbind(combDat, curData)
  
  }
  colnames(combDat) = c("Sample_id", "MetaseqCount", "AmpseqCount")
  return(combDat)
}

sumDat = sumTable(df = ludy_df, minDepth = 10)

sum(sumDat$MetaseqCount) + sum(sumDat$AmpseqCount)

write.csv(sumDat, "SeqSamplesNumbers_Depth10.csv", row.names = F)