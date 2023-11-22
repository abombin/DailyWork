

library(viridis)
library(stringr)
ampseq = read.csv("Paper/Tables/Original_Multivar.csv")

ampseq = ampseq[!ampseq$Predictor == "Intercept",]

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

ampseq<-addPCat(dataTab=ampseq, pcol = "Pval", newCol = "P.value")
# ampseq$Full_variable_name = str_to_title(ampseq$Full_variable_name)
# ampseq$Variable = toupper(ampseq$Variable)

p1 <- ggplot(ampseq, aes(x = Estimate, y  = reorder(Predictor, Estimate), fill = Pval)) +
  geom_bar(stat = "identity", color = "black") +
  labs( x = "Estimate", y = "Predictor") + theme_classic() + theme(text = element_text(size = 26), plot.title = element_text(hjust = 0.5, vjust = 0.5))

ggsave(file = "Paper/Figs/Shannon_Multivar.png", plot = p1, height = 14, width = 18, units = "in", dpi = 300)