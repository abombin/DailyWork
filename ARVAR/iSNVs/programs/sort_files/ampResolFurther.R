
missAmp = read.csv("amp_uniqueMiss.csv")

# can be parallel
statFailedLib = function(inDir) {
  dirsList = list.files(inDir)
  outDir = "amp_miss_stats/"
  dir.create(outDir, recursive = T, showWarnings = F)
  for ( curDir in dirsList) {
    cmd_str = paste0("seqkit stats -a ", inDir, "/", curDir)
    curStats = system(cmd_str, intern = T)
    outFile = paste0(outDir, curDir)
    write.table(curStats, outFile, row.names = F, col.names = F, quote = F)
  }
}

statFailedLib(inDir = "amp_miss_files")

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

ampStatsMiss = combStats(inDir = "amp_miss_stats")
colnames(ampStatsMiss)[17] = "fileName"
missAmp$fileName = basename(missAmp$Pbnas_path)

combMiss = plyr::join(missAmp, ampStatsMiss, by = "fileName", type = "left", match = "all")

combMissSel = combMiss[, c("combSeqName", "combOrigMiss", "Pbnas_path", "File_size" , "Spike.in.confirmed", "min_len", "avg_len", "max_len")]

colnames(combMissSel)[1:2] = c("Real_Lib_Name", "Original_Lib_Name")

write.csv(combMissSel, "amp_miss_resolved_andrei.csv", row.names = F)
system("aws s3 cp amp_miss_resolved_andrei.csv s3://abombin/ARVAR/iSNVs/")