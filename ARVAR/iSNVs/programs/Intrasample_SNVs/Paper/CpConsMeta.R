library(doParallel)
library(foreach)


metaseq = read.csv("Paper/metaseqWu_derep_decont_v3.csv")

samplesList = unique(metaseq$OrigName)

getFilesList = function(inPath, pattern) {
  filesList = list.files(inPath, pattern = pattern, full.names = T)
  return(filesList)
}

metaList = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/output_2/variants/bowtie2", pattern = ".sorted.bam$")
meta1 = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/output/variants/bowtie2", pattern = ".sorted.bam$")
meta2 = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-29/output/variants/bowtie2", pattern = ".sorted.bam$")
#meta3 = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/EHC-C19-2120P_S30_L001/output_EHC-C19-2120P_S30_L001/variants/bowtie2", pattern = ".sorted.bam$")
meta_comb = c(metaList, meta1, meta2)

samplesPaths = character()

for ( i in samplesList) {
  cursample = paste0(i, ".sorted.bam")
  curPath = meta_comb[grepl(cursample, meta_comb)]
  samplesPaths = c(samplesPaths, curPath)
}

length(unique(samplesPaths))

baseMeta = basename(samplesPaths)
freqDf = data.frame(table(baseMeta))
baseMeta = gsub(".sorted.bam", "",baseMeta )
baseMeta = sub(".*EHC", "EHC", baseMeta)
length(unique(baseMeta))
samplesList[!samplesList%in%baseMeta]
freqDf = data.frame(table(baseMeta))

samplesPaths[grepl("EHC-C19-1636Z_S16_L001", samplesPaths)]
samplesList[grepl("EHC-C19-1636Z_S16_L001", samplesList)]

samplesPaths <- samplesPaths[samplesPaths != "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/output_2/variants/bowtie2/p21175-s016_EHC-C19-1636Z_S16_L001.sorted.bam"]
length(unique(samplesPaths))


consPath = gsub("bowtie2", "ivar/consensus/bcftools", samplesPaths)
consPath = gsub(".sorted.bam", ".consensus.fa", consPath)


dir.create("Consens_fasta/Metaseq/", recursive = T)

prepRefs = function(curPath) {
  sampleName = basename(curPath)
  sampleName = gsub(".consensus.fa", "", sampleName)
  outName = paste0("Consens_fasta/Metaseq/", sampleName, ".fa")
  file.copy(curPath, outName)
  
}

cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = consPath) %dopar% {
  prepRefs(curPath = i)
}
parallel::stopCluster(cl = cl)