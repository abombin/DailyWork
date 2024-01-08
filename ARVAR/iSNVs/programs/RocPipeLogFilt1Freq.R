library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(plyr)


getConsensus <- function(metaSeq, ampSeq, protocol, maxFreq, minFreq, freqCol) {
  metaseq <- read.csv(metaSeq)
  ampseq <- read.csv(ampSeq)
  metaseq =  metaseq[metaseq$Coverage >= 97,]
  ampseq =   ampseq[ampseq$Coverage >= 97,]
  
  metaseqFilt <- metaseq[metaseq[[freqCol]] >= minFreq & metaseq[[freqCol]] <= maxFreq, ]
  ampseqFilt <- ampseq[ampseq[[freqCol]] >= minFreq & ampseq[[freqCol]] <= maxFreq, ]
  
  ampseqFilt = ampseqFilt[ampseqFilt$Sample%in%metaseqFilt$Sample,]
  metaseqFilt = metaseqFilt[metaseqFilt$Sample%in%ampseqFilt$Sample,]
  
  
  if (protocol == "ampseq") {
    # dataframe with stats
    targDf = ampseqFilt
    
    # data to check agains
    refDf = metaseqFilt
    targSnv <- unique(refDf$Samp_Pos_Ref_Alt)
  } else if (protocol == "metaseq") {
    # dataframe with stats
    targDf = metaseqFilt
    
    # data to check agains
    refDf = ampseqFilt
    targSnv <- unique(refDf$Samp_Pos_Ref_Alt)
  }
  
  ConsTest <- numeric()
  
  for (i in 1:nrow(targDf)) {
    curSnv =  targDf[i, "Samp_Pos_Ref_Alt"]
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


# check filtering 
# q1 = ampseqFreqCons[, c("ALLELE.FREQUENCY", "refFreq" )]
# q2 = q1[q1$ALLELE.FREQUENCY < 0.02 & q1$refFreq < 0.02, ]
# q3 = q1[q1$ALLELE.FREQUENCY > 0.98 & q1$refFreq > 0.98, ]

runRoc = function(df, protocol, freqCol, splitPerc) {
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  
  set.seed(42)
  train_idx <- createDataPartition(dfFilt$ConsTest, p = splitPerc, list = FALSE)
  train_data <- dfFilt[train_idx, ]
  test_data <- dfFilt[-train_idx, ]
  
  if (protocol == "metaseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + I(Mean_depth^2)" ))
    aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    
  } else if (protocol == "ampseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + I(Mean_depth^2)"))
    aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
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
  return(test_data)
}

ampseqFreqCons = getConsensus(metaSeq="snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_0_v2.csv", ampSeq="snvs_comb_res/ampseq_overlap_comb_derep_covFilt_0_v2.csv", 
                               protocol="ampseq", maxFreq=0.98, minFreq=0.02, freqCol="Freq_adj")

#ampseqFreqCons = ampseqFreqCons[order(ampseqFreqCons$Samp_Pos_Ref_Alt),]

df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")
df = df[df$ALLELE.FREQUENCY <= 0.98,]
# df = df[order(df$Samp_Pos_Ref_Alt),]
# identical(df$Samp_Pos_Ref_Alt,ampseqFreqCons$Samp_Pos_Ref_Alt )

ampseqFreqConsPred = runRoc(df=ampseqFreqCons, protocol = "ampseq", freqCol = "Freq_adj", splitPerc = 0.7)
table(ampseqFreqCons$ConsTest)


# metaseq
metaseqFreqCons = getConsensus(metaSeq="snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_97.csv", ampSeq="snvs_comb_res/ampseq_overlap_comb_derep_covFilt_97.csv", 
                              protocol="metaseq", maxFreq=0.98, minFreq=0.02, freqCol="ALLELE.FREQUENCY")

df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")
df = df[df$ALLELE.FREQUENCY <= 0.98,]

metaseqFreqConsPred = runRoc(df=metaseqFreqCons , protocol = "metaseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)
table(metaseqFreqCons$ConsTest)
table(df$ConsTest)

