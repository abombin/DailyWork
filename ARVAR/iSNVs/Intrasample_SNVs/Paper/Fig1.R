library(ggplot2)

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




