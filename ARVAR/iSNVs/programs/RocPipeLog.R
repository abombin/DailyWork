library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(plyr)


filterFreq = function(targDf, refDf, maxFreq, minFreq, freqCol) {
  refDfSel = refDf[, c("Samp_Pos_Ref_Alt", "ALLELE.FREQUENCY", "Freq_adj")]
  colnames(refDfSel)[2:3] = c("refFreq", "refFreqAdj")
  
  targDfCons = targDf[targDf$ConsTest == 1,]
  targDfSpur = targDf[targDf$ConsTest == 0,]
  
  targDfCons = plyr::join(targDfCons, refDfSel, by = "Samp_Pos_Ref_Alt", type = "left", match = "all")
  if (freqCol == "ALLELE.FREQUENCY") {
    targDfConsFilt = targDfCons[targDfCons[[freqCol]] <= maxFreq | targDfCons$refFreq <= maxFreq , ]
    targDfConsFilt = targDfConsFilt[targDfConsFilt[[freqCol]] >= minFreq | targDfConsFilt$refFreq  >= minFreq ,]
    
  } else if (freqCol == "Freq_adj") {
    targDfConsFilt = targDfCons[targDfCons[[freqCol]] <= maxFreq | targDfCons$refFreqAdj <= maxFreq , ]
    targDfConsFilt = targDfConsFilt[targDfConsFilt[[freqCol]] >= minFreq | targDfConsFilt$refFreqAdj  >= minFreq ,]
  }
  # filter spurious calls by frequency
  targDfSpur =  targDfSpur[ targDfSpur[[freqCol]] >= minFreq &  targDfSpur[[freqCol]] <= maxFreq, ]
  #targDfSpur =  targDfSpur[ targDfSpur[[freqCol]] >= minFreq,]
  
  dfFilt = rbind.fill(targDfConsFilt, targDfSpur)
  
  return(dfFilt)
}


getConsensus <- function(metaSeq, ampSeq, protocol, maxFreq, minFreq, freqCol) {
  metaseq <- read.csv(metaSeq)
  ampseq <- read.csv(ampSeq)
  
  #metaseqFilt <- metaseq[metaseq[[freqCol]] >= minFreq & metaseq[[freqCol]] <= maxFreq, ]
  #ampseqFilt <- ampseq[ampseq[[freqCol]] >= minFreq & ampseq[[freqCol]] <= maxFreq, ]
  
  if (protocol == "ampseq") {
    # dataframe with stats
    targDf = ampseq
    
    # data to check agains
    refDf = metaseq
    targSnv <- refDf$Samp_Pos_Ref_Alt
  } else if (protocol == "metaseq") {
    # dataframe with stats
    targDf = metaseq
    
    # data to check agains
    refDf = ampseq
    targSnv <- refDf$Samp_Pos_Ref_Alt
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
  
  filtDf = filterFreq(targDf=targDf, refDf=refDf, maxFreq=maxFreq, minFreq=minFreq, freqCol=freqCol)
  
  return(filtDf)
}


# check filtering 
# q1 = ampseqFreqCons[, c("ALLELE.FREQUENCY", "refFreq" )]
# q2 = q1[q1$ALLELE.FREQUENCY < 0.02 & q1$refFreq < 0.02, ]
# q3 = q1[q1$ALLELE.FREQUENCY > 0.98 & q1$refFreq > 0.98, ]

runRoc = function(dfFilt, protocol, freqCol, splitPerc) {
  dfFilt = dfFilt[!is.na(dfFilt$Var_Al_RelPos),]
  
  set.seed(42)
  train_idx <- createDataPartition(dfFilt$ConsTest, p = splitPerc, list = FALSE)
  train_data <- dfFilt[train_idx, ]
  test_data <- dfFilt[-train_idx, ]
  
  if (protocol == "metaseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + I(Mean_depth^2)" ))
    aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    
  } else if (protocol == "ampseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos + I(Mean_depth^2)"))
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

ampseqFreqCons = getConsensus(metaSeq="snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_97.csv", ampSeq="snvs_comb_res/ampseq_overlap_comb_derep_all.csv", 
                              protocol="ampseq", maxFreq=0.98, minFreq=0.02, freqCol="ALLELE.FREQUENCY")

ampseqFreqConsPred = runRoc(dfFilt=ampseqFreqCons, protocol = "ampseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)
table(ampseqFreqCons$ConsTest)




metaseqFreqCons = getConsensus(metaSeq="snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_97.csv", ampSeq="snvs_comb_res/ampseq_overlap_comb_derep_covFilt_97.csv", 
                              protocol="metaseq", maxFreq=1, minFreq=0.02, freqCol="ALLELE.FREQUENCY")

ampseqFreqConsPred = runRoc(dfFilt=metaseqFreqCons, protocol = "metaseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)
table(metaseqFreqCons$ConsTest)
