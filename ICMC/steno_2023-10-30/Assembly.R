library(foreach)
library(doParallel)

filesList = list.files("processPar")

trimGalore = function(curFile) {
  setwd(paste0("processPar/", curFile))
  r1 = paste0(curFile, "_R1_001.fastq.gz")
  r2 = paste0(curFile, "_R2_001.fastq.gz")
  cmd_str = paste0("/home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20 -j 2 --basename sample  --paired ", r1, " ", r2)
  system(cmd_str)
  setwd("../../")
}

cl = makeCluster(4, type = "FORK")
registerDoParallel(cl)

combRes = foreach(i = filesList) %dopar% {
  #trimGalore(curFile = i)
}

stopCluster(cl)

dir.create("gtdbtk/input/", recursive = T, showWarnings = F)


Spades = function(curFile) {
  setwd(paste0("processPar/", curFile))
  cmd_str = paste0("/home/ubuntu/spades/SPAdes-3.15.5-Linux/bin/spades.py -1 sample_val_1.fq.gz -2 sample_val_2.fq.gz -o spades_out -t 8 --isolate")
  #system(cmd_str)
  
  system(paste0("cp spades_out/contigs.fasta ../../gtdbtk/input/", curFile, ".fasta"))
  
  setwd("../../")
}

cl = makeCluster(4, type = "FORK")
registerDoParallel(cl)

combRes = foreach(i = filesList) %dopar% {
  Spades(curFile = i)
}

stopCluster(cl)