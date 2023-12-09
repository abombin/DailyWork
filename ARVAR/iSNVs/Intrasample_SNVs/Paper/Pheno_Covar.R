# Format phenotype data
# get dist and mds covars


library(phylotools)

snvsDat = read.csv("Paper/Additional_analysis/GWAS/AmpMetaCombSnvsPostPredDedup_Wu.csv")

fasta = read.fasta("Paper/allAlignment/allCons.align")

ampseqLength = length(list.files('Consens_fasta/Ampseq/'))
metaseqLength = length(list.files('Consens_fasta/Metaseq/'))
ampseqLength + metaseqLength + 1 == nrow(fasta)

protAmpseq = rep("Ampseq", ampseqLength)
protMetaseq = rep("Metaseq", metaseqLength)
protocol = c("Wuhan", protAmpseq, protMetaseq)

fasta$Protocol = protocol

fasta$seq.name = gsub("\\ .*", "", fasta$seq.name)

# filter alignment