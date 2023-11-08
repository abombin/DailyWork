library(ggplot2)

# 
# df = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv')
# df = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv')
# 
# sumTable = data.frame(table(df$All_Alignment_Pos, df$VAR.NT))
# colnames(sumTable)[1:2] = c("Position", "Variant")
# 
# sumTable = sumTable[sumTable$Freq > 10,]
# 
# length(unique(sumTable$Position))
# 
# 
# 
# p1 =  ggplot(sumTable, aes(x = Position, y = Freq, label = Variant)) +
#   geom_text(size = 3, vjust = -0.5) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme(text = element_text(size = 28))
# 
# 
# ggsave(file = "metaseq_freq10.png", plot=p1, height = 20, width = 40, dpi = 300, units = "in", limitsize = FALSE)



df1 = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv')
df2 = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv')

df = rbind(df1, df2)

sumTable = data.frame(table(df$All_Alignment_Pos, df$VAR.NT))
colnames(sumTable)[1:2] = c("Position", "Variant")

sumTable = sumTable[sumTable$Freq > 50,]

length(unique(sumTable$Position))



p1 =  ggplot(sumTable, aes(x = Position, y = Freq, label = Variant)) +
  geom_text(size = 3, vjust = -0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(text = element_text(size = 28))


ggsave(file = "combine_freq50.png", plot=p1, height = 20, width = 40, dpi = 300, units = "in", limitsize = FALSE)


length(unique(df1$AllAlignPos_Var))
length(unique(df2$AllAlignPos_Var))

length(unique(df1$AllAlignPos_Var[df1$AllAlignPos_Var %in% df2$AllAlignPos_Var]))
