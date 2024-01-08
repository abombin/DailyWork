failAmp = read.csv("amp_uniqueFailFiltr.csv")

# can be parallel
statFailedLib = function(inDir) {
  dirsList = list.files(inDir)
  outDir = "amp_fail_stats/"
  dir.create(outDir, recursive = T, showWarnings = F)
  for ( curDir in dirsList) {
    cmd_str = paste0("seqkit stats -a ", inDir, "/", curDir)
    curStats = system(cmd_str, intern = T)
    outFile = paste0(outDir, curDir)
    write.table(curStats, outFile, row.names = F, col.names = F, quote = F)
  }
}

combStats = function(inDir) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  filesList = list.files(inDir)
  for ( curFile in filesList ) {
    #print(curFile)
    inFile = paste0(inDir, "/", curFile)
    curDf = read.table(inFile, T)
    curDf$NewDir = curFile
    combDat = rbind(combDat, curDf)
  }
  return(combDat)
}

statFailedLib(inDir = "amp_fail_files")
ampStatsFail = combStats(inDir = "amp_fail_stats")
colnames(ampStatsFail)[17] = "fileName"
failAmp$fileName = basename(failAmp$Pbnas_path)

combFail = plyr::join(failAmp, ampStatsFail, by = "fileName", type = "left", match = "all")

combFailSel = combFail[, c("combSeqName", "Pbnas_path", "File_size" , "Spike.in.confirmed", "min_len", "avg_len", "max_len")]

colnames(combFailSel)[1] = c("Real_Lib_Name")

write.csv(combFailSel, "amp_fail_resolved_andrei.csv", row.names = F)
system("aws s3 cp amp_fail_resolved_andrei.csv s3://abombin/ARVAR/iSNVs/")

identical(combFailSel$min_len, combFailSel$max_len)