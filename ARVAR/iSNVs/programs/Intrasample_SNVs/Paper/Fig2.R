
library(viridis)
library(stringr)
ampseq = read.csv("Paper/Tables/Ampseq_Multivar_LogisticCoef_Table.csv")

ampseq = ampseq[!ampseq$Full_variable_name == "Intercept",]

addPCat<-function(dataTab, pcol, newCol){
  dataTab[newCol]<-NA
  for (i in 1:nrow(dataTab)){
    if (dataTab[[pcol]][i] > 0.05){
      dataTab[[newCol]][i]<- " > 0.05"
    } else if (dataTab[[pcol]][i] < 0.05 & dataTab[[pcol]][i] >= 0.01  ) {
      dataTab[[newCol]][i]<- " < 0.05"
    } else if (dataTab[[pcol]][i] < 0.01 & dataTab[[pcol]][i] >= 0.001) {
      dataTab[[newCol]][i]<- " < 0.01"
    } else if (dataTab[[pcol]][i] < 0.001) {
      dataTab[[newCol]][i]<- " < 0.001"
    }
  }
  return(dataTab)
}

ampseq<-addPCat(dataTab=ampseq, pcol = "Pvalue", newCol = "P.value")
ampseq$Full_variable_name = str_to_title(ampseq$Full_variable_name)
ampseq$Variable = toupper(ampseq$Variable)

p1 <- ggplot(ampseq, aes(x = Estimate, y  = reorder(Variable, Estimate), fill = P.value)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Amplicon Sequencing", x = "Estimate", y = "Variables") + theme_classic() + theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5, vjust = 0.5))

ggsave(file = "Paper/Figs/Ampseq_Multivar_LogisticCoef.png", plot = p1, height = 14, width = 18, units = "in", dpi = 300)

# metaseq
metaseq = read.csv("Paper/Tables/Metaseq_Multivar_LogisticCoef_Table.csv")
metaseq = metaseq[!metaseq$Full_variable_name == "Intercept",]
metaseq<-addPCat(dataTab=metaseq, pcol = "Pvalue", newCol = "P.value")
metaseq$Full_variable_name = str_to_title(metaseq$Full_variable_name)
metaseq$Variable = toupper(metaseq$Variable)

p2 <- ggplot(metaseq, aes(x = Estimate, y  = reorder(Variable, Estimate), fill = P.value)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Metagenome Sequencing", x = "Estimate", y = "Variables") + theme_classic() + theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5, vjust = 0.5))

ggsave(file = "Paper/Figs/Metaseq_Multivar_LogisticCoef.png", plot = p2, height = 14, width = 18, units = "in", dpi = 300)

# random forest
ampseq = read.csv("IntraSnv_results/SumModels/Ampseq_VarImportance.csv")
ampseq$variable = toupper(ampseq$variable)
p3 <- ggplot(ampseq, aes(x = importance, y  = reorder(variable, importance))) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Amplicon Sequencing", x = "Importance", y = "Variables") + theme_classic() + theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5, vjust = 0.5))

ggsave(file = "Paper/Figs/Ampseq_RandForImportance.png", plot = p3, height = 14, width = 18, units = "in", dpi = 300)

metaseq = read.csv("IntraSnv_results/SumModels/Metaseq_VarImportance.csv")
metaseq$variable = toupper(metaseq$variable)
p4 <- ggplot(metaseq, aes(x = importance, y  = reorder(variable, importance))) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Metagenome Sequencing", x = "Importance", y = "Variables") + theme_classic() + theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5, vjust = 0.5))

ggsave(file = "Paper/Figs/Metaseq_RandForImportance.png", plot = p4, height = 14, width = 18, units = "in", dpi = 300)
