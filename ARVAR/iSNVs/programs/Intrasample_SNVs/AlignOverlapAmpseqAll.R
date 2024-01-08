# metaseq processing
library(doParallel)
library(foreach)

metaseq = read.csv("snvs_comb_res/ampseq_all_v2_names.csv")
samplesList = unique(metaseq$OrigName)
length(unique(samplesList))


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

prepFastqs = function(curPath) {
  sampleName = basename(curPath)
  sampleName = gsub(".ivar_trim.sorted.bam", "", sampleName)
  targDir = paste0("IntraSnv_ampseq_all/", sampleName, "/")
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

cl <- makeCluster(16, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = samplesPaths) %dopar% {
  prepFastqs(curPath = i)
}
parallel::stopCluster(cl = cl)

consPath = gsub("bowtie2", "ivar/consensus/bcftools", samplesPaths)
consPath = gsub(".ivar_trim.sorted.bam", ".consensus.fa", consPath)


prepRefs = function(curPath) {
  sampleName = basename(curPath)
  sampleName = gsub(".consensus.fa", "", sampleName)
  targDir = paste0("IntraSnv_ampseq_all/", sampleName, "/")
  outFasta = paste0(targDir, "reference.fa")
  file.copy(curPath, outFasta)
  ref_cmd = paste0("hisat2-build -p 4 ", outFasta, " ", targDir,  "reference")
  system(ref_cmd)
  align_str = paste0("hisat2 -p 4 -x ", targDir, "reference -U ", targDir, "output.fastq -S ", targDir, "output.sam")
  system(align_str)
  file.remove("output.fastq")
  bam_str = paste0("samtools view -bu -F 4 ", targDir, "output.sam -o - | samtools sort - -o ", targDir, "output.bam")
  #system(bam_str)
  ind_str = paste0("samtools index ", targDir, "output.bam")
  #system(ind_str)
}

cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = consPath) %dopar% {
  prepRefs(curPath = i)
}
parallel::stopCluster(cl = cl)