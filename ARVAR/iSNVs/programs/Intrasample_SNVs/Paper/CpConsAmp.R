library(doParallel)
library(foreach)


ampseq = read.csv("Paper/ampseqWu_derep_v3.csv")

samplesList = unique(ampseq$OrigName)

getFilesList = function(inPath, pattern) {
  filesList = list.files(inPath, pattern = pattern, full.names = T)
  return(filesList)
}

ampList = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_ampseq_2023-08-22/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
amp1 = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
amp2 = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
amp_comb = c(ampList, amp1, amp2)

samplesPaths = character()

for ( i in samplesList) {
  cursample = paste0(i, ".ivar_trim.sorted.bam")
  curPath = amp_comb[grepl(cursample, amp_comb)]
  samplesPaths = c(samplesPaths, curPath)
}

length(unique(samplesPaths))

baseMeta = basename(samplesPaths)
freqDf = data.frame(table(baseMeta))

samplesPaths[grepl("EHC-C19-1444P-LAmp_S35_L001.ivar_trim.sorted.bam", samplesPaths)]
samplesList[grepl("EHC-C19-1444P-LAmp_S35_L001", samplesList)]

samplesPaths <- samplesPaths[samplesPaths != "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/output/variants/bowtie2/EHC-C19-1444P-LAmp_S35_L001.ivar_trim.sorted.bam"]
length(unique(samplesPaths))

consPath = gsub("bowtie2", "ivar/consensus/bcftools", samplesPaths)
consPath = gsub(".ivar_trim.sorted.bam", ".consensus.fa", consPath)

dir.create("Consens_fasta/Ampseq/", recursive = T)

prepRefs = function(curPath) {
  sampleName = basename(curPath)
  sampleName = gsub(".consensus.fa", "", sampleName)
  outName = paste0("Consens_fasta/Ampseq/", sampleName, ".fa")
  file.copy(curPath, outName)

}

cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = consPath) %dopar% {
  prepRefs(curPath = i)
}
parallel::stopCluster(cl = cl)