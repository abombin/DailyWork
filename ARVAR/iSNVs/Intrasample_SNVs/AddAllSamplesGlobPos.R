adjustPositions = function(df, protocol, targDir) {
  combDat = data.frame()
  sampleNames = unique(df$OrigName) 
  for (curName in sampleNames) {
    dfSub = df[df$OrigName == curName,]
    if (protocol=="ampseq") {
      indPath = paste0(targDir, "/ampseq/", curName, ".csv")
    } else if (protocol=="metaseq") {
      indPath = paste0(targDir, "/metaseq/", curName, ".csv")
    }
    try({
      indDf = read.csv(indPath)
      indDf = indDf[, c("Original_Pos", "Alignment_Pos")]
      colnames(indDf)[1:2] = c("POSITION", "All_Alignment_Pos")
      joinDat = plyr::join(dfSub, indDf, by = "POSITION", type = "left", match = "all")
      joinDat$AllAlignPos_Var = paste(joinDat$All_Alignment_Pos, joinDat$VAR.NT, sep = "__")
      joinDat$AllAlignPos_Ref_Var = paste(joinDat$All_Alignment_Pos, joinDat$REF.NT,  joinDat$VAR.NT,  sep = "__")
      joinDat$Sample_AllAlignPos_Ref_Var = paste(joinDat$Sample, joinDat$All_Alignment_Pos, joinDat$REF.NT, joinDat$VAR.NT, sep = "__")
      combDat = rbind(combDat, joinDat)
    })
  }
  return(combDat)
}

# ampseq
ampseq = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv')
min(ampseq$coverage)
ampseqPos = adjustPositions(df=ampseq, protocol="ampseq", targDir = "All_Samples_Mapping/AllSampAlign_Indecies/")
sumTable = data.frame(table(ampseqPos$AllAlignPos_Ref_Var))
unique(is.na(ampseqPos$AllAlignPos_Ref_Var))
write.csv(ampseqPos, "IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv", row.names = F)

#metaseq
metaseq = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions.csv')
min(metaseq$coverage)
metaseqPos = adjustPositions(df=metaseq, protocol="metaseq", targDir = "All_Samples_Mapping/AllSampAlign_Indecies/")
sumTableMeta = data.frame(table(metaseqPos$AllAlignPos_Ref_Var))
unique(is.na(metaseqPos$AllAlignPos_Ref_Var))
write.csv(metaseqPos, "IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv", row.names = F)

system("aws s3 cp IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv s3://abombin/ARVAR/iSNVs/IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv")

system("aws s3 cp IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv s3://abombin/ARVAR/iSNVs/IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv")