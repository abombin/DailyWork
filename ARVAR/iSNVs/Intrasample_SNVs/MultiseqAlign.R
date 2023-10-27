

dir.create("All_Samples_Mapping/metaseq/Comb_fasta/", recursive = T, showWarnings = F)
dir.create("All_Samples_Mapping/ampseq/Comb_fasta/", recursive = T, showWarnings = F)

makeAlignment = function(combStats, protocol) {
  system(paste0("cat references/MN908947.3.fna > All_Samples_Mapping/", protocol, "/Comb_fasta/allSequences.fasta"))
  for ( i in 1:nrow(combStats) ) {
    ampSample = combStats$OrigName[i]
    metaSample = combStats$OrigName[i]
    sampName = combStats$Sample[i]
    catFile = paste0("All_Samples_Mapping/", protocol, "/Comb_fasta/allSequences.fasta")
    ampPath = paste0("IntraSnv_ampseq_all/",  ampSample, "/reference.fa")
    metaPath = paste0("IntraSnv_metaseq_all/",  metaSample, "/reference.fa")
    if (protocol == "ampseq" ) {
      cmd_str = paste("cat ", ampPath, ">>", catFile, sep = " ")
    } else if (protocol == "metaseq" ) {
      cmd_str = paste("cat ", metaPath, ">>", catFile, sep = " ")
    }
    system(cmd_str)
  }
  cmd_mafft = paste0("mafft --maxiterate 1000 --thread 70 --globalpair ",  catFile, " > All_Samples_Mapping/", protocol, "/All_samples_align.fasta")
  system(cmd_mafft)
  
}

# run alignment function

# metaseq = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions.csv')
# metaseqSamp = data.frame(unique(metaseq[, c("Sample", "OrigName")]))
# 
# makeAlignment(combStats=metaseqSamp, protocol = 'metaseq')


ampseq = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv')
ampseqSamp = data.frame(unique(ampseq[, c("Sample", "OrigName")]))

makeAlignment(combStats=ampseqSamp, protocol = 'ampseq')