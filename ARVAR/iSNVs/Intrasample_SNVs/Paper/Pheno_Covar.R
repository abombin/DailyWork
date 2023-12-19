# Format phenotype data
# get dist and mds covars

setwd("/home/ubuntu/extraVol/ARVAR/iSNVs")

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
fasta$OrigName_Prot = paste(fasta$seq.name, fasta$Protocol, sep = "__")

fastaFilt = fasta[fasta$OrigName_Prot%in%snvsDat$OrigName_Prot,]
nrow(fastaFilt) == length(unique(snvsDat$Sample))
length(unique(fastaFilt$seq.name)) == nrow(fastaFilt)

selMeta = unique(snvsDat[, c("Sample", "OrigName", "OrigName_Prot")])

combFasta = plyr::join(fastaFilt, selMeta, by= "OrigName_Prot", type = "left", match = "all")
colnames(combFasta)

combFasta = combFasta[, c("Sample", "seq.text")]
colnames(combFasta)[1] = "seq.name" 

setwd("/home/ubuntu/extraVol/ARVAR/iSNVs/Paper/Additional_analysis/GWAS/")

dat2fasta(combFasta, "selAlignGwasWu.fasta")

system("../../../snp-dists/snp-dists selAlignGwasWu.fasta > selAlignGwasWu.dist")

distMat =  read.delim("selAlignGwasWu.dist")

colnames(distMat)[1] = "Sample"
rownames(distMat) = distMat$Sample
distMat = distMat[, 2:ncol(distMat)]
colnames(distMat) = rownames(distMat)

mds_result <- cmdscale(distMat, eig = TRUE, k = 10)

# Extract eigenvalues
eigenvalues <- mds_result$eig

# Create a scree plot
scree_plot <- plot(1:10, eigenvalues[1:10], type = "b", 
                   xlab = "Number of Components", ylab = "Eigenvalues",
                   main = "Scree Plot")

mds_components <- data.frame(mds_result$points[, 1:5])
mds_components$Sample = rownames(mds_components)
write.csv(mds_components, "Mds_covar_Wu.csv", row.names = F)


# try with principal components
pca_result <- prcomp(distMat, center = TRUE, scale. = TRUE)

prop_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Check the length of prop_var
cat("Length of prop_var:", length(prop_var), "\n")

# Plot the scree plot
scree_plot <- plot(1:10, prop_var[1:10],
                   type = "b", xlab = "Principal Component Number", ylab = "Proportion of Variance",
                   main = "Scree Plot for the First 10 Principal Components")
