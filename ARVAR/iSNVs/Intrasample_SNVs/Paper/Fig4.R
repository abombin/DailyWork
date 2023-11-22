ampseq = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv')

sumTab = data.frame(table(ampseq$All_Alignment_Pos))

sumTab$Position = as.numeric(as.character(sumTab$Var1))

p1 =  ggplot(sumTab, aes(x = Position, y = Freq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 28))+
  theme(text = element_text(size = 28)) +
  scale_x_continuous(breaks = seq(0, max(sumTab$Position), by = 2000))

ggsave(file = "Paper/Figs/FreqPlot_SNVs_Ampseq.png", plot = p1, height = 12, width = 18, units = "in", dpi = 300)

# metaseq
metaseq = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv')
sumTab = data.frame(table(metaseq$All_Alignment_Pos))
sumTab$Position = as.numeric(as.character(sumTab$Var1))

p2 =  ggplot(sumTab, aes(x = Position, y = Freq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 28))+
  theme(text = element_text(size = 28)) +
  scale_x_continuous(breaks = seq(0, max(sumTab$Position), by = 2000))

ggsave(file = "Paper/Figs/FreqPlot_SNVs_Metaseq.png", plot = p2, height = 12, width = 18, units = "in", dpi = 300)

combDat = rbind(ampseq, metaseq)

sumTab = data.frame(table(combDat$All_Alignment_Pos))
sumTab$Position = as.numeric(as.character(sumTab$Var1))

p3 =  ggplot(sumTab, aes(x = Position, y = Freq)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 28))+
  theme(text = element_text(size = 28)) +
  scale_x_continuous(breaks = seq(0, max(sumTab$Position), by = 2000))

ggsave(file = "Paper/Figs/FreqPlot_SNVs_Comb.png", plot = p3, height = 12, width = 18, units = "in", dpi = 300)