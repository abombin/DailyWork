library(ggvenn)
library(ggpubr)
library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)
library(ggplot2)

ampseq = read.csv("IntraSnv_results/ampseq_comb_derep_overlap.csv")
metaseq = read.csv("IntraSnv_results/metaseq_comb_derep_overlap.csv")

ampseq = ampseq[ampseq$coverage >= 97,]
metaseq = metaseq[metaseq$coverage >= 97,]

ampseqOverlap = ampseq[ampseq$Sample%in%metaseq$Sample,]
metaseqOverlap = metaseq[metaseq$Sample %in% ampseq$Sample,]

length(unique(ampseqOverlap$OrigName))
length(unique(metaseqOverlap$OrigName))

adjustPositions = function(df, protocol, targDir) {
  combDat = data.frame()
  sampleNames = unique(df$OrigName) 
  for (curName in sampleNames) {
    dfSub = df[df$OrigName == curName,]
    if (protocol=="ampseq") {
      indPath = paste0(targDir, "/Ampseq/", curName, ".csv")
    } else if (protocol=="metaseq") {
      indPath = paste0(targDir, "/Metaseq/", curName, ".csv")
    }
    try({
      indDf = read.csv(indPath)
      indDf = indDf[, c("Original_Pos", "Alignment_Pos")]
      colnames(indDf)[1] = "POSITION"
      joinDat = plyr::join(dfSub, indDf, by = "POSITION", type = "left", match = "all")
      joinDat$Sample_AlignPos_Var = paste(joinDat$Sample, joinDat$Alignment_Pos, joinDat$VAR.NT, sep = "__")
      joinDat$AlignPos_Var = paste(joinDat$Alignment_Pos, joinDat$VAR.NT, sep = "__")
      joinDat$Sample_AlignPos_Ref_Var = paste(joinDat$Sample, joinDat$Alignment_Pos, joinDat$REF.NT, joinDat$VAR.NT, sep = "__")
      combDat = rbind(combDat, joinDat)
    })
  }
  return(combDat)
}

getConsensus <- function(metaSeq, ampSeq, protocol, maxFreq, minFreq, freqCol, Snv_col) {
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
    targSnv <- unique(refDf[, Snv_col])
  } else if (protocol == "metaseq") {
    # dataframe with stats
    targDf = metaseqFilt
    # data to check agains
    refDf = ampseqFilt
    targSnv <- unique(refDf[, Snv_col])
  }
  ConsTest <- numeric()
  for (i in 1:nrow(targDf)) {
    curSnv =  targDf[i, Snv_col]
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
    aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 500)
  } else if (protocol == "ampseq") {
    # glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + I(meandepth^2)"))
    # aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + meandepth"))
    glm_formula = as.formula("ConsTest ~ ALLELE.FREQUENCY+ STRAND.BIAS + QUAL + Var_Al_RelPos + Ref_Al_RelPos + meandepth + coverage + meanmapq + meanbaseq")
    aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 500)
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
    ggtitle(curTitle)+
    theme(text = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 24))
  ggsave(filename =  paste0("Venn/Intra_SampAlignPosRefVar_AllAlign_2_", curTitle,'.jpeg'), plot = venPlot, width = 9, height = 9, units = 'in', dpi = 600, device = 'jpeg')
}




ampseqOverlap = adjustPositions(df=ampseqOverlap, protocol="ampseq", targDir = "Overlap_Pos_mapping_2/AllSampAlign_Indecies")
metaseqOverlap = adjustPositions(df=metaseqOverlap, protocol="metaseq", targDir = "Overlap_Pos_mapping_2/AllSampAlign_Indecies") 

ampseqCons = getConsensus(metaSeq=metaseqOverlap, ampSeq=ampseqOverlap, protocol="ampseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY", Snv_col = "Sample_AlignPos_Ref_Var")
metaseqCons = getConsensus(metaSeq=metaseqOverlap, ampSeq=ampseqOverlap, protocol="metaseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY", Snv_col = "Sample_AlignPos_Ref_Var")

# write.csv(ampseqCons, "IntraSnv_results/ampseq_ConsTest_freq_1_0.01.csv", row.names = F)
# write.csv(metaseqCons, "IntraSnv_results/metaseq_ConsTest_freq_1_0.01.csv", row.names = F)

metaseqFreq1 = metaseqCons[metaseqCons$ALLELE.FREQUENCY == 1,]
ampseqFreq1 = ampseqCons[ampseqCons$ALLELE.FREQUENCY == 1,]


table(ampseqCons$ConsTest)
table(metaseqCons$ConsTest)

getVenn(metaseq=metaseqCons, ampseq=ampseqCons, filtSteps=1, maxFreq = 1, minFreq = 0, Snv_col="Sample_AlignPos_Ref_Var")

runRoc(df=metaseqCons, protocol="metaseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)
runRoc(df=ampseqCons, protocol="ampseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)

length(unique(ampseqCons$Sample_AlignPos_Ref_Var))
length(unique(metaseqCons$Sample_AlignPos_Ref_Var))

q1 = data.frame(table(metaseqCons$Sample_AlignPos_Ref_Var))
q2 = metaseqCons[metaseqCons$Sample_AlignPos_Var == "EHC-C19-1424V__8921__C",]