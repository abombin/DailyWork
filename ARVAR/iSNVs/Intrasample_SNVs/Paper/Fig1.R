library(ggplot2)
library(ggvenn)

ampseq = read.csv('IntraSnv_results/ampseq_comb_derep.csv')
metaseq = read.csv('IntraSnv_results/metaseq_comb_derep.csv')

ampseq = ampseq[ampseq$coverage >= 97,]
metaseq = metaseq[metaseq$coverage >= 97, ]

ampseq$Protocol = "Ampseq"
metaseq$Protocol = "Metaseq"

combDat = rbind(ampseq, metaseq)

sumDat = data.frame(table(combDat$Protocol, combDat$OrigName))
sumDat = sumDat[sumDat$Freq > 0 , ]

colnames(sumDat) = c("Protocol", "Sample", "Frequency")

p1 <- ggplot(sumDat, aes(x=Protocol, y=Frequency, fill = Protocol)) + 
  geom_violin(trim=FALSE)+
  geom_jitter() + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 8,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

ggsave(file = "Paper/Figs/AllSamples_SNVs_Violin.png", plot = p1, height = 10, width = 14, units = "in", dpi = 300)


# frequency of SNVs
result <- aggregate(ALLELE.FREQUENCY ~ Protocol + OrigName, data = combDat, FUN = mean)


ggplot(result, aes(x=Protocol, y=ALLELE.FREQUENCY, fill = Protocol)) + 
  geom_boxplot() + theme_classic() + theme(text = element_text(size = 26))

p2 <- ggplot(result, aes(x=Protocol, y=ALLELE.FREQUENCY, fill = Protocol)) + 
  geom_violin(trim=FALSE)+
  geom_jitter() + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 8,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

ggsave(file = "Paper/Figs/AllSamples_SNVs_MeanFreq_Violin.png", plot = p2, height = 10, width = 14, units = "in", dpi = 300)

# make Venn diagrams

getVenn<-function(metaseq, ampseq, filtSteps, maxFreq, minFreq, Snv_col){
  ampSnps<-ampseq[, Snv_col]
  # process meta
  metaSnps<-metaseq[, Snv_col]
  print(length(metaSnps))
  print(length(ampSnps))
  print(length(metaSnps[metaSnps%in%ampSnps]))
  print(length(metaSnps[metaSnps%in%ampSnps]) / (length(metaSnps) + length(ampSnps)) * 100 )
  # make a list
  snps<-list('Amplicon'=ampSnps, 'Metagenomic'=metaSnps)
  curTitle = paste0("FiltSteps_", filtSteps, "_MaxFreq_", maxFreq, "_minFreq_", minFreq)
  # plot
  venPlot<-ggvenn(snps, text_size=5, set_name_size=8)+
    #ggtitle(curTitle)+
    theme(text = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 24))
  ggsave(filename =  paste0("Paper/Figs/Venn_OverlapSamples_RandFor.png"), plot = venPlot, width = 9, height = 9, units = 'in', dpi = 600, device = 'jpeg')
}

ampseq = read.csv("IntraSnv_results/ampseq_ConsTest_freq_1_0.csv")
metaseq = read.csv("IntraSnv_results/metaseq_ConsTest_freq_1_0.csv")

getVenn(metaseq=metaseq, ampseq=ampseq, filtSteps=1, maxFreq = 1, minFreq = 0, Snv_col="Sample_AlignPos_Ref_Var")
