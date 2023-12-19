ampseq = read.csv('IntraSnv_results/ampseq_comb_derep.csv')
metaseq = read.csv('IntraSnv_results/metaseq_comb_derep.csv')

ampseq = ampseq[ampseq$coverage >= 97,]
metaseq = metaseq[metaseq$coverage >= 97, ]

ampseqSamples = unique(ampseq$OrigName)
metaseqSamples = unique(metaseq$OrigName)

ampPaths = paste0("IntraSnv_ampseq_all/", ampseqSamples, "/reference.fa")
metaPaths = paste0("IntraSnv_metaseq_all/", metaseqSamples, "/reference.fa")


dir.create("Paper/allAlignment/Comb_fasta/", recursive = T, showWarnings = F)

system("cat references/MN908947.3.fna > Paper/allAlignment/Comb_fasta/allSequences.fasta")

for (i in ampPaths) {
  system(paste0("cat ", i, " >> Paper/allAlignment/Comb_fasta/allSequences.fasta"))
}

for (i in metaPaths) {
  system(paste0("cat ", i, " >> Paper/allAlignment/Comb_fasta/allSequences.fasta"))
}

system("mafft --maxiterate 1000 --thread 70 --globalpair Paper/allAlignment/Comb_fasta/allSequences.fasta > Paper/allAlignment/allSamples.align")