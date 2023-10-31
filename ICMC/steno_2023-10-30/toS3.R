makeCustomOut = function() {
    newDir = paste0("custom_output/")
    dir.create(newDir, recursive = T, showWarnings = F)
    curDirs = list.files(paste0("panOut/bactopia-runs"), full.names = T)
    # change next line when run with normal bactopia runs
    for (curDir in curDirs) {
      distFile = paste0(curDir, "/core-genome.distance.tsv")
      treeFile = paste0(curDir, "/core-genome.aln.gz.treefile")
      contTree = paste0(curDir, "/core-genome.aln.gz.contree")
      iqFile = paste0(curDir, "/core-genome.aln.gz.iqtree")
      try({
        system(paste0("cp ", distFile, " ", newDir))
        system(paste0("cp ", treeFile, " ", newDir))
        system(paste0("cp ", contTree, " ", newDir))
        system(paste0("cp ", iqFile, " ", newDir))
      })
      
    }
}

getS3Path = function() {
  curPath = getwd()
  curList = strsplit(curPath, "/")
  dirPath = paste(c(curList[[1]][5], curList[[1]][6], curList[[1]][7]), collapse = "/")
  s3Path = paste(c("s3://transfer-files-emory", dirPath, "custom_output/"), collapse = "/" )
  return(s3Path)
}


transferS3 = function(curDir, s3Path) {
  cmd_str = paste0("aws s3 cp --recursive ", curDir, '/ ', s3Path)
  system(cmd_str)
}

makeCustomOut()
s3Path = getS3Path()
transferS3(curDir = "custom_output", s3Path = s3Path)

getS3Path = function() {
  curPath = getwd()
  curList = strsplit(curPath, "/")
  dirPath = paste(c(curList[[1]][5], curList[[1]][6], curList[[1]][7]), collapse = "/")
  s3Path = paste(c("s3://transfer-files-emory", dirPath), collapse = "/" )
  return(s3Path)
}

s3Path = getS3Path()

bactopiaOut = paste0(s3Path, "/", "bactopia_output/")
panOut = paste0(paste0(s3Path, "/", "pangenome/"))

cmd_str = paste0("aws s3 cp --recursive bactopia_output/ ", bactopiaOut)
system(cmd_str)

cmd_str = paste0("aws s3 cp --recursive panOut/ ", panOut)
system(cmd_str)