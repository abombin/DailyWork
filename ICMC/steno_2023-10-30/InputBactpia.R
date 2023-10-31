getType_001 = function(inPath, runType) {
  filesList = list.files(inPath, pattern = "_R1_001.fastq.gz")
  if (length(filesList) > 0) {
    sample = gsub('_R1_001.fastq.gz', "", filesList )
    r1 = normalizePath(paste0(inPath, "/",  filesList))
    r2 = normalizePath(paste0(inPath, "/", sample, "_R2_001.fastq.gz"))
    runtype = runType
    genome_size = 0
    species = "Stenotrophomonas maltophilia"
    df = data.frame(sample, runtype, genome_size, species, r1, r2)
    df$extra = ""
    
  } 
  else {
    df = NULL
  }
  return(df)
}


prevSamples = getType_001(inPath = "prev_run_fastqs", runType = "paired-end")

curSamples = getType_001(inPath = "fastqs", runType = "paired-end")


combInput = rbind(prevSamples, curSamples)

write.table(combInput, "bactopiaInput.tsv", row.names = F, col.names = T, quote = F, sep = "\t")