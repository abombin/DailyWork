library(ggvenn)
library(ggpubr)
library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)
library(ggplot2)

ampseq = read.csv("IntraSnv_ampseq_overlap/ampseq_comb_stats.csv")
metaseq = read.csv("IntraSnv_metaseq_overlap/metaseq_comb_stats.csv")

ampseq = ampseq[ampseq$coverage >= 97,]
metaseq = metaseq[metaseq$coverage >= 97,]

ampseqFilt = ampseq[ampseq$Sample%in%metaseq$Sample,]
metaseqFilt = metaseq[metaseq$Sample %in% ampseq$Sample,]

length(unique(ampseqFilt$OrigName))
length(unique(metaseqFilt$OrigName))

adjustPositions = function(df, protocol) {
  combDat = data.frame()
  sampleNames = unique(df$OrigName) 
  for (curName in sampleNames) {
    dfSub = df[df$OrigName == curName,]
    if (protocol=="ampseq") {
      indPath = paste0("Overlap_Pos_mapping/Indecies/Ampseq/", curName, ".csv")
    } else if (protocol=="metaseq") {
      indPath = paste0("Overlap_Pos_mapping/Indecies/Metaseq/", curName, ".csv")
    }
    try({
      indDf = read.csv(indPath)
      indDf = indDf[, c("Original_Pos", "Alignment_Pos")]
      colnames(indDf)[1] = "POSITION"
      joinDat = plyr::join(dfSub, indDf, by = "POSITION", type = "left", match = "all")
      joinDat$Sample_AlignPos_Var = paste(joinDat$Sample, joinDat$Alignment_Pos, joinDat$VAR.NT, sep = "__")
      #joinDat$AlignPos_Var = paste(joinDat$Alignment_Pos, joinDat$VAR.NT, sep = "__")
      combDat = rbind(combDat, joinDat)
    })
  }
  return(combDat)
}

getConsensus <- function(metaSeq, ampSeq, protocol, maxFreq, minFreq, freqCol) {
    metaseq <- metaSeq
    ampseq <- ampSeq
    metaseq =  metaseq[metaseq$coverage >= 97,]
    ampseq =   ampseq[ampseq$coverage >= 97,]
    metaseqFilt <- metaseq[metaseq[[freqCol]] >= minFreq & metaseq[[freqCol]] <= maxFreq, ]
    ampseqFilt <- ampseq[ampseq[[freqCol]] >= minFreq & ampseq[[freqCol]] <= maxFreq, ]
    ampseqFilt = ampseqFilt[ampseqFilt$Sample%in%metaseqFilt$Sample,]
    metaseqFilt = metaseqFilt[metaseqFilt$Sample%in%ampseqFilt$Sample,]
    if (protocol == "ampseq") {
      # dataframe with stats
      targDf = ampseqFilt
      # data to check agains
      refDf = metaseqFilt
      targSnv <- unique(refDf$Sample_AlignPos_Var)
    } else if (protocol == "metaseq") {
      # dataframe with stats
      targDf = metaseqFilt
      # data to check agains
      refDf = ampseqFilt
      targSnv <- unique(refDf$Sample_AlignPos_Var)
    }
    ConsTest <- numeric()
    for (i in 1:nrow(targDf)) {
      curSnv =  targDf[i, "Sample_AlignPos_Var"]
      if (curSnv%in%targSnv){
        curCons = 1
      } else {
        curCons = 0
      }
      ConsTest = c(ConsTest, curCons)
    }
    targDf$ConsTest <- ConsTest
    return(targDf)
}

runRoc = function(df, protocol, freqCol, splitPerc) {
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  set.seed(42)
  train_idx <- createDataPartition(dfFilt$ConsTest, p = splitPerc, list = FALSE)
  train_data <- dfFilt[train_idx, ]
  test_data <- dfFilt[-train_idx, ]
  if (protocol == "metaseq") {
    # glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + I(meandepth^2)" ))
    # aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    glm_formula <- as.formula(paste("ConsTest ~ ", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos  + meandepth" ))
    glm_formula = as.formula("ConsTest ~ ALLELE.FREQUENCY+ STRAND.BIAS + QUAL + Var_Al_RelPos + Ref_Al_RelPos + meandepth + coverage + meanmapq + meanbaseq")
    aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 700)
  } else if (protocol == "ampseq") {
    # glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + I(meandepth^2)"))
    # aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + meandepth"))
    glm_formula = as.formula("ConsTest ~ ALLELE.FREQUENCY+ STRAND.BIAS + QUAL + Var_Al_RelPos + Ref_Al_RelPos + meandepth + coverage + meanmapq + meanbaseq")
    aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 700)
  }
  probs <- predict(aucModel, newdata = test_data, type = "response")
  roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
  # actual predictions and testing
  threshold <- 0.5
  test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)
  equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)
  # Count the number of times values in Column1 are NOT equal to Column2
  not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)
  sumCorrect = equal_count / (equal_count + not_equal_count) * 100
  print(roc_obj$auc)
  print(sumCorrect )
}


getVenn<-function(metaseq, ampseq, filtSteps, maxFreq, minFreq){
  ampSnps<-ampseq$Sample_AlignPos_Var
  # process meta
  metaSnps<-metaseq$Sample_AlignPos_Var
  print(length(metaSnps))
  print(length(ampSnps))
  print(length(metaSnps[metaSnps%in%ampSnps]))
  print(length(metaSnps[metaSnps%in%ampSnps]) / (length(metaSnps) + length(ampSnps)) * 100 )
  # make a list
  snps<-list('Amplicon'=ampSnps, 'Metagenomic'=metaSnps)
  curTitle = paste0("FiltSteps_", filtSteps, "_MaxFreq_", maxFreq, "_minFreq_", minFreq)
  # plot
  venPlot<-ggvenn(snps, text_size=5, set_name_size=8)+
    ggtitle(curTitle)+
    theme(text = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 24))
  #ggsave(filename =  paste0("Venn/Intra_AlignPos_", curTitle,'.jpeg'), plot = venPlot, width = 9, height = 9, units = 'in', dpi = 600, device = 'jpeg')
}




ampseqFilt = adjustPositions(df=ampseqFilt, protocol="ampseq")
metaseqFilt = adjustPositions(df=metaseqFilt, protocol="metaseq") 

ampseqCons = getConsensus(metaSeq=metaseqFilt, ampSeq=ampseqFilt, protocol="ampseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY")
metaseqCons = getConsensus(metaSeq=metaseqFilt, ampSeq=ampseqFilt, protocol="metaseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY")

getVenn(metaseq=metaseqCons, ampseq=ampseqCons, filtSteps=1, maxFreq = 1, minFreq = 0)

runRoc(df=metaseqCons, protocol="metaseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)
runRoc(df=ampseqCons, protocol="ampseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)


# write.csv(metaseqCons, "IntraSnv_metaseq_overlap/metaseq_ampseq_overlap_97_allFreq.csv", row.names= F)
# write.csv(ampseqCons, "IntraSnv_ampseq_overlap/ampseq_metaseq_overlap_97_allFreq.csv", row.names = F)

# make histogram ggplot
ggplot(ampseqCons, aes(x=ALLELE.FREQUENCY)) + geom_histogram() + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  ggtitle("Ampseq Freq 0.01-1") + theme(text = element_text(size = 26))
ggplot(metaseqCons, aes(x=ALLELE.FREQUENCY)) + geom_histogram() + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) + ggtitle("Metaseq Freq 0.01-1") +  theme(text = element_text(size = 26))

# try without indels

ampseqCons = getConsensus(metaSeq=metaseqFilt, ampSeq=ampseqFilt, protocol="ampseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY")
metaseqCons = getConsensus(metaSeq=metaseqFilt, ampSeq=ampseqFilt, protocol="metaseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY")

getVenn<-function(metaseq, ampseq, filtSteps, maxFreq, minFreq){
  ampSnps<-ampseq$Sample_AlignPos_Var
  # process meta
  metaSnps<-metaseq$Sample_AlignPos_Var
  print(length(metaSnps))
  print(length(ampSnps))
  print(length(metaSnps[metaSnps%in%ampSnps]))
  print(length(metaSnps[metaSnps%in%ampSnps]) / (length(metaSnps) + length(ampSnps)) * 100 )
  # make a list
  snps<-list('Amplicon'=ampSnps, 'Metagenomic'=metaSnps)
  curTitle = paste0("FiltSteps_", filtSteps, "_MaxFreq_", maxFreq, "_minFreq_", minFreq)
  # plot
  venPlot<-ggvenn(snps, text_size=5, set_name_size=8)+
    ggtitle(curTitle)+
    theme(text = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 24))
  ggsave(filename =  paste0("Venn/Intra_AlignPos_ExclIndels_", curTitle,'.jpeg'), plot = venPlot, width = 9, height = 9, units = 'in', dpi = 600, device = 'jpeg')
}

metaseqSub = metaseqCons[nchar(metaseqCons$REF.NT) == 1 & nchar(metaseqCons$VAR.NT) == 1, ]
ampseqSub = ampseqCons[nchar(ampseqCons$REF.NT) == 1 & nchar(ampseqCons$VAR.NT) == 1, ]

runRoc(df=metaseqSub, protocol="metaseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)
runRoc(df=ampseqSub, protocol="ampseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)

getVenn(metaseq=metaseqSub, ampseq=ampseqSub, filtSteps=1, maxFreq = 1, minFreq = 0)



