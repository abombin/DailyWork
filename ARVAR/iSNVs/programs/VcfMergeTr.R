dir.create('Paper/Additional_analysis/GWAS/VcfWu/', recursive = T)

ampseq = read.csv("Paper/ampseqWu_derep_v3.csv")

ampseqSamp = unique(ampseq[, c("OrigName", "Batch")])

# this way seems to make all samples look like one sample
for (i in 1:nrow(ampseqSamp)) {
  curBatch = ampseqSamp$Batch[i]
  curSamp = ampseqSamp$OrigName[i]
  if (curBatch == "New") {
     vcfPath = paste0("ampseq_vivacity_found/", curSamp, "/sample_lf.vcf")
     system(paste0("bgzip ", vcfPath))
     vcfPath = paste0("ampseq_vivacity_found/", curSamp, "/sample_lf.vcf.gz")
     system(paste0("tabix -p vcf ",  vcfPath))
     outPath = paste0("Paper/Additional_analysis/GWAS/VcfWu/", curSamp, "_ampseq.vcf.gz")
     system(paste0("bcftools view ", vcfPath, " -O z -o ", outPath))
     system(paste0("tabix -p vcf ",  outPath))
  }
}


for (i in 1:nrow(ampseqSamp)) {
  curBatch = ampseqSamp$Batch[i]
  curSamp = ampseqSamp$OrigName[i]
  if (curBatch == "New") {
    vcfPath = paste0("ampseq_vivacity_found/", curSamp, "/sample_lf.vcf")
    outPath = paste0("Paper/Additional_analysis/GWAS/VcfWu/", curSamp, "_ampseq.vcf")
    file.copy(vcfPath, outPath)
    system(paste0("bgzip ", outPath))
    
    system(paste0("tabix -p vcf ",  outPath, ".gz"))
  }
}


for (i in 1:nrow(ampseqSamp)) {
  curBatch = ampseqSamp$Batch[i]
  curSamp = ampseqSamp$OrigName[i]
  if (curBatch == "New") {
    vcfPath = paste0("ampseq_vivacity_found/", curSamp, "/sample_lf.vcf")
    outPath = paste0("Paper/Additional_analysis/GWAS/VcfWu/", curSamp, "_ampseq.vcf.gz")
    exp_str = paste0("bcftools norm -O z -m - ", vcfPath, " > ", outPath)
    system(exp_str)
    
    system(paste0("tabix -p vcf ", outPath))
  }
}


rewriteVcf = function(input_file, output_file, curSamp) {
  con_in <- file(input_file, "r")
  con_out <- file(output_file, "w")
  
  # Read and process each line
  while (length(line <- readLines(con_in, n = 1)) > 0) {
    if (startsWith(line, "#")) {
      # If line starts with #, write it as is to the output file
      writeLines(line, con_out)
    } else {
      # Replace "MN908947.3" with "MySample1" and write to the output file
      modified_line <- gsub("MN908947\\.3", curSamp, line)
      writeLines(modified_line, con_out)
    }
  }
  
  # Close connections
  close(con_in)
  close(con_out)
}

rewriteVcf = function(input_file, output_file, curSamp) {
  con_in <- file(input_file, "r")
  con_out <- file(output_file, "w")
  
  # Read and process each line
  while (length(line <- readLines(con_in, n = 1)) > 0) {
      # Replace "MN908947.3" with "MySample1" and write to the output file
      modified_line <- gsub("MN908947\\.3", curSamp, line)
      writeLines(modified_line, con_out)

  }
  
  # Close connections
  close(con_in)
  close(con_out)
}

for (i in 1:nrow(ampseqSamp)) {
  curBatch = ampseqSamp$Batch[i]
  curSamp = ampseqSamp$OrigName[i]
  if (curBatch == "New") {
    vcfPath = paste0("ampseq_vivacity_found/", curSamp, "/sample_lf.vcf")
    outPath = paste0("Paper/Additional_analysis/GWAS/VcfWu/", curSamp, "_ampseq.vcf")
    
    rewriteVcf(input_file=vcfPath, output_file=outPath, curSamp=curSamp)
    
    system(paste0("bgzip ", outPath))
    
    system(paste0("tabix -p vcf ",  outPath, ".gz"))
  }
}



# seems that vcf files has to be compressed with bgzip and indexed
system("bcftools merge -m none -0 -O v Paper/Additional_analysis/GWAS/VcfWu/*.vcf.gz > Paper/Additional_analysis/GWAS/mergedSNVWu.vcf")