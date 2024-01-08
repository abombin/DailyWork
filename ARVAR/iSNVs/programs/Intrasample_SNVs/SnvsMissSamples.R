library(doParallel)
library(foreach)

listDirs = list.files("IntraSnv_ampseq_all")

listOutput = list.files("IntraSnv_ampseq_all", pattern = "sample_lofreq-output.tsv", recursive = T)

getMissing = function(listDirs, listOutput) {
  missingDirs = character()
  for ( i in listDirs ) {
    if (!any(grepl(i, listOutput)) ) {
      missingDirs = c(missingDirs, i)
    }
  }
  return(missingDirs)
}

missingOutput = getMissing(listDirs, listOutput)

