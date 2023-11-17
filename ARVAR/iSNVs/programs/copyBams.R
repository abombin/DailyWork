library(doParallel)
library(foreach)

# copy newly found ampseq bam files
filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_ampseq_2023-08-22/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
samplesList = gsub(".ivar_trim.sorted.bam", "", filesList)

copyBams = function(curSample, inDir, outDir, protocol) {
  curDir = paste0(outDir, "/", curSample, "/")
  dir.create(curDir, showWarnings = F, recursive = T)
  if (protocol == "ampseq") {
    curFiles = list.files(inDir, pattern = paste0(curSample, ".*ivar_trim.sorted.bam"), full.names = T)
  } else {
    curFiles = list.files(inDir, pattern = curSample, full.names = T)
  }
  for (curFile in curFiles) {
    if (grepl(".bai", curFile)) {
      file.copy(curFile, paste0(curDir, "output.bam.bai"))
    } else {
      file.copy(curFile, paste0(curDir, "output.bam"))
    }
  }
}

runPar = function(samplesList, inDir, outDir, protocol) {
  cores = 93
  curCluster = makeCluster(cores, type = "FORK")
  registerDoParallel(curCluster)
  curRun = foreach(i = samplesList) %dopar% {
    copyBams(curSample = i, inDir = inDir, 
             outDir = outDir, protocol = protocol)
  }
  parallel::stopCluster(curCluster)
}

#runPar(samplesList = samplesList, inDir = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_ampseq_2023-08-22/output/variants/bowtie2", outDir = "ampseq_vivacity_found", protocol = "ampseq")

# copy ampseq previously found samples
# batch 1 
filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
samplesList = gsub(".ivar_trim.sorted.bam", "", filesList)
#runPar(samplesList = samplesList, inDir = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/output/variants/bowtie2", outDir = "ampseq_vivacity_old", protocol = "ampseq")
# batch 2
filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
samplesList = gsub(".ivar_trim.sorted.bam", "", filesList)
#runPar(samplesList = samplesList, inDir = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/output/variants/bowtie2", outDir = "ampseq_vivacity_old", protocol = "ampseq")

##  Copy newly found metaseq

filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/output/variants/bowtie2", pattern = ".sorted.bam$")
samplesList = gsub(".sorted.bam", "", filesList)
runPar(samplesList = samplesList, inDir = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/output/variants/bowtie2", outDir = "metaseq_vivacity_found", protocol = "metaseq")

# copy previously found metaseq
# batch1
dir.create("metaseq_vivacity_old", showWarnings = F)
filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/output/variants/bowtie2", pattern = ".sorted.bam$")
samplesList = gsub(".sorted.bam", "", filesList)
runPar(samplesList = samplesList, inDir = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/output/variants/bowtie2", outDir = "metaseq_vivacity_old", protocol = "metaseq")
# batch2
filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-29/output/variants/bowtie2", pattern = ".sorted.bam$")
samplesList = gsub(".sorted.bam", "", filesList)
runPar(samplesList = samplesList, inDir = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-29/output/variants/bowtie2", outDir = "metaseq_vivacity_old", protocol = "metaseq")