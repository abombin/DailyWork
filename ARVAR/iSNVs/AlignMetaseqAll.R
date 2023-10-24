# metaseq processing
library(doParallel)
library(foreach)

metaseq = read.csv("snvs_comb_res/metaseq_all_v2_names.csv")
samplesList = unique(metaseq$OrigName)
length(unique(samplesList))


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

samplesPaths <- samplesPaths[samplesPaths != "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/output/variants/bowtie2/p21175-s016_EHC-C19-1636Z_S16_L001.sorted.bam"]
length(unique(samplesPaths))

prepFastqs = function(curPath) {
  sampleName = basename(curPath)
  print(curPath)
  sampleName = gsub(".sorted.bam", "", sampleName)
  targDir = paste0("IntraSnv_metaseq_all/", sampleName, "/")
  dir.create(targDir, recursive = T, showWarnings = F)
  #cmd_str = paste0("samtools view -b -F 4 ", curPath, " | hisat2 -x references/MN908947.3 -U - -S - | samtools view -b -o ", targDir, "output.bam -")
  #cmd_str = paste0("samtools view -h -F 4 ", curPath, " > ", targDir, "aligned_reads.sam")
  cmd_str = paste0("samtools view -b -F 4 ", curPath, " | bedtools bamtofastq -i - -fq ", targDir, "output.fastq")
  system(cmd_str)
  # align_str = paste0("hisat2 -x references/MN908947.3 -U ", targDir, "output.fastq -S ", targDir, "output.sam")
  # system(align_str)
  # bam_str = paste0("samtools view -bu -F 4 ", targDir, "output.sam -o - | samtools sort - -o ", targDir, "output.bam")
  # system(bam_str)
  
}

cl <- makeCluster(30, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = samplesPaths) %dopar% {
  prepFastqs(curPath = i)
}
parallel::stopCluster(cl = cl)

consPath = gsub("bowtie2", "ivar/consensus/bcftools", samplesPaths)
consPath = gsub(".sorted.bam", ".consensus.fa", consPath)


prepRefs = function(curPath) {
  sampleName = basename(curPath)
  sampleName = gsub(".consensus.fa", "", sampleName)
  targDir = paste0("IntraSnv_metaseq_all/", sampleName, "/")
  outFasta = paste0(targDir, "reference.fa")
  file.copy(curPath, outFasta)
  ref_cmd = paste0("hisat2-build -p 4 ", outFasta, " ", targDir,  "reference")
  system(ref_cmd)
  align_str = paste0("hisat2 -p 4 -x ", targDir, "reference -U ", targDir, "output.fastq -S ", targDir, "output.sam")
  system(align_str)
  bam_str = paste0("samtools view -bu -F 4 ", targDir, "output.sam -o - | samtools sort - -o ", targDir, "output.bam")
  #system(bam_str)
  ind_str = paste0("samtools index ", targDir, "output.bam")
  #system(ind_str)
  file.remove("output.fastq")
}

cl <- makeCluster(8, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = consPath) %dopar% {
  prepRefs(curPath = i)
}
parallel::stopCluster(cl = cl)