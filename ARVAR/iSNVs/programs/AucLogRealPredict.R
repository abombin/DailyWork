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

makeRealPredictions = function(trainData, testData, freqCol) {
  set.seed(42)
  train_data = trainData[!is.na(trainData$Var_Al_RelPos),]
  test_data = testData[!is.na(testData$Var_Al_RelPos),]
  test_data = test_data[!test_data$Samp_Pos_Ref_Alt%in%train_data$Samp_Pos_Ref_Alt,]
  glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + I(Mean_depth^2)" ))
  aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
  # review stats on the same data
  probs <- predict(aucModel, newdata = train_data, type = "response")
  roc_obj <- roc(train_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
  threshold <- 0.5
  train_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)
  equal_count <- sum(train_data$ConsTest == train_data$PredictedConsTest)
  not_equal_count <- sum(train_data$ConsTest != train_data$PredictedConsTest)
  sumCorrect = equal_count / (equal_count + not_equal_count) * 100
  print(roc_obj$auc)
  print(sumCorrect )
  # make real predictions
  threshold <- 0.5
  probs <- predict(aucModel, newdata = test_data, type = "response")
  test_data$ConsTest <- ifelse(probs >= threshold, 1, 0)
  return(test_data)

}

combineData = function(realDat, predictDat) {
  # Convert ConsTest to numeric
  realDat$ConsTest = as.numeric(realDat$ConsTest)
  predictDat$ConsTest = as.numeric(predictDat$ConsTest)
  # Filter realDat where ConsTest is 1
  realDatFilt = realDat[realDat$ConsTest == 1, ]
  # Remove rows with NA in realDatFilt$ConsTest
  realDatFilt = realDatFilt[!is.na(realDatFilt$ConsTest), ]
  # Filter predictDat where ConsTest is 1
  predictDatFilt = predictDat[predictDat$ConsTest == 1, ]
  # Remove rows with NA in predictDatFilt$ConsTest
  predictDatFilt = predictDatFilt[!is.na(predictDatFilt$ConsTest), ]
  # Combine filtered data frames
  combDat = rbind(realDatFilt, predictDatFilt)
  # Print unique values of ConsTest
  print(unique(combDat$ConsTest))
  return(combDat)
}

cols_to_scale <- c("ALLELE.FREQUENCY", "QUAL", "STRAND.BIAS", "Var_Al_RelPos", "Mean_depth")

df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97_v2.csv")
#df1 = df[df$ALLELE.FREQUENCY <= 0.98,]
dfFilt = filterFreq(df=df, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)
dfFilt$Origina_Mean_Depth = dfFilt$Mean_depth
dfFilt$Original_freq = dfFilt$ALLELE.FREQUENCY
# scale
dfFilt[, cols_to_scale] <- scale(dfFilt[, cols_to_scale])
metaseqFreqConsPred = runRoc(df=dfFilt , protocol = "metaseq", freqCol = "ALLELE.FREQUENCY", splitPerc = 0.7)

all_snvs = read.csv("snvs_comb_res/metaseq_comb_derep_decont_covFilt_97_v2.csv")
all_snvs_filt = filterFreq(df=all_snvs, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)
all_snvs_filt$Origina_Mean_Depth = all_snvs_filt$Mean_depth
all_snvs_filt$Original_freq = all_snvs_filt$ALLELE.FREQUENCY
# Scale the selected columns and assign them back to the data frame
all_snvs_filt[, cols_to_scale] <- scale(all_snvs_filt[, cols_to_scale])

predictedDat = makeRealPredictions(trainData=dfFilt, testData=all_snvs_filt, freqCol="ALLELE.FREQUENCY")

#table(predictedDat$ConsTest)
# 90571  1511 
# table(dfFilt$ConsTest)
# 9504  953 

combDat = combineData(realDat=dfFilt, predictDat=predictedDat)
table(predictedDat$combDat)

write.csv(combDat, "snvs_comb_res/metaseq_LogPredict_scale_covFilt_97_v2.csv", row.names = F)
