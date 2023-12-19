# process overlapping samples
metaseq = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")

metaseqOverlapSamples = metaseq[metaseq$ConsTest == 1,]
metaseqOverlapSamples = data.frame(unique(metaseqOverlapSamples$Sample))
colnames(metaseqOverlapSamples) = "Sample"

# process metaseq samples
metaseq = read.csv("snvs_comb_res/metaseq_comb_derep_decont_covFilt_97.csv")
metaseqSamples = unique(metaseq[,c("Sample", "OrigName")])

#process ampseq samples
ampseq = read.csv("snvs_comb_res/ampseq_comb_derep_covFilt_97.csv")
ampseqSamples = unique(ampseq[,c("Sample", "OrigName")])

##
write.csv(metaseqOverlapSamples, "model_summary/overlap_samples_97_coverage.csv", row.names = F)
write.csv(metaseqSamples, "model_summary/metaseq_samples_97_coverage.csv", row.names = F)
write.csv(ampseqSamples, "model_summary/ampseq_samples_97_coverage.csv", row.names = F)