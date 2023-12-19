
library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)


filterFreq = function(df, freqCol, maxFreq, minFreq) {
  ampseq = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")
  metaseq = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")
  
  ampseqFilt = ampseq[ampseq[[freqCol]] >= minFreq & ampseq[[freqCol]] <= maxFreq, ]
  print(nrow(ampseqFilt) == length(unique(ampseqFilt$Samp_Pos_Ref_Alt)))
  metaseqFilt = metaseq[metaseq[[freqCol]] >= minFreq & metaseq[[freqCol]] <= maxFreq, ]
  print(nrow(metaseqFilt) == length(unique(metaseqFilt$Samp_Pos_Ref_Alt)))
  selSnps = unique(c(ampseqFilt$Samp_Pos_Ref_Alt, metaseqFilt$Samp_Pos_Ref_Alt))
  
  dfFilt = df[df$Samp_Pos_Ref_Alt%in%selSnps,]
  return(dfFilt)
}

runRoc = function(df, protocol, freqCol, splitPerc) {
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  
  set.seed(42)
  train_idx <- createDataPartition(dfFilt$ConsTest, p = splitPerc, list = FALSE)
  train_data <- dfFilt[train_idx, ]
  test_data <- dfFilt[-train_idx, ]
  
  if (protocol == "metaseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + Mean_depth" ))
    aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 3000)
    
  } else if (protocol == "ampseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + Mean_depth"))
    aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 3000)
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

#metaseq
df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")
dfFilt = filterFreq(df=df, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)
metaseqFreqConsPred = runRoc(df=dfFilt , protocol = "metaseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)

dfFilt = filterFreq(df=df, freqCol="Freq_adj", maxFreq=0.98, minFreq=0.02)
metaseqFreqConsPred = runRoc(df=dfFilt , protocol = "metaseq", freqCol = "Freq_adj", splitPerc = 0.7)


## ampseq
df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")
dfFilt = filterFreq(df=df, freqCol="Freq_adj", maxFreq=0.98, minFreq=0.02)
ampseqFreqConsPred = runRoc(df=dfFilt , protocol = "ampseq", freqCol = "Freq_adj", splitPerc = 0.7)

dfFilt = filterFreq(df=df, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)
ampseqFreqConsPred = runRoc(df=dfFilt , protocol = "ampseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)