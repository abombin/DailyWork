getFilesList = function(inPath) {
  filesList = list.files(inPath)
  sampleNames = gsub("\\_R1_.*", "", filesList)
  sampleNames = gsub("\\_R2_.*", "", sampleNames)
  sampleNames = unique(sampleNames)
  return(sampleNames)
}

filesList = getFilesList(inPath = "fastqs")

copyFiles = function(filesList, inPath) {
  for ( i in filesList ) {
    curFiles = list.files(inPath, full.names = T, pattern = i)
    targDir = paste0("processPar/", i, "/")
    dir.create(targDir, recursive = T, showWarnings = T)
    for (curFile in curFiles) {
      system(paste0("cp ", curFile, " ", targDir))
    }
  }
}

copyFiles(filesList=filesList, inPath=inPath)

