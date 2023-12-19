
library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)


filterFreq = function(df, freqCol, maxFreq, minFreq) {
  dfFilt = df[df[[freqCol]] >= minFreq & df[[freqCol]] <= maxFreq, ]
  return(dfFilt)
}

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



df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97_v2.csv")
#df1 = df[df$ALLELE.FREQUENCY <= 0.98,]
dfFilt = filterFreq(df=df, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)
metaseqFreqConsPred = runRoc(df=dfFilt , protocol = "metaseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)

dfFilt = filterFreq(df=df, freqCol="Freq_adj", maxFreq=0.98, minFreq=0.02)
metaseqFreqConsPred = runRoc(df=dfFilt , protocol = "metaseq", freqCol = "Freq_adj", splitPerc = 0.7)


## ampseq
df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97_v2.csv")
dfFilt = filterFreq(df=df, freqCol="Freq_adj", maxFreq=0.98, minFreq=0.02)
ampseqFreqConsPred = runRoc(df=dfFilt , protocol = "ampseq", freqCol = "Freq_adj", splitPerc = 0.7)

dfFilt = filterFreq(df=df, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)
ampseqFreqConsPred = runRoc(df=dfFilt , protocol = "ampseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)

