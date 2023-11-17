
library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)


# Freq_adj  ALLELE.FREQUENCY

filterFreq = function(df, freqCol, maxFreq, minFreq) {
  dfFilt = df[df[[freqCol]] >= minFreq & df[[freqCol]] <= maxFreq, ]
  return(dfFilt)
}

df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")
dfFilt = filterFreq(df=df, freqCol="ALLELE.FREQUENCY", maxFreq=0.98, minFreq=0.02)

set.seed(42)

dfFilt = dfFilt[!is.na(dfFilt$Var_Al_RelPos),]
#dfFilt$ConsTest = as.factor(dfFilt$ConsTest)

train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]


rfModel <- randomForest(ConsTest  ~ ALLELE.FREQUENCY   + QUAL + Var_Al_RelPos   + Mean_depth , data = train_data , ntree = 1500)
# Predict probabilities
probs <- predict(rfModel, newdata = test_data, type = "response")
#probs = as.numeric(probs)
# Calculate ROC
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
threshold <- 0.5
test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)
equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)
# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)
sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print(sumCorrect )

## Freq_adj


df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")
dfFilt = filterFreq(df=df, freqCol="Freq_adj", maxFreq=0.98, minFreq=0.02)

set.seed(42)

dfFilt = dfFilt[!is.na(dfFilt$Var_Al_RelPos),]
#dfFilt$ConsTest = as.factor(dfFilt$ConsTest)

train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]


rfModel <- randomForest(ConsTest  ~ Freq_adj   + QUAL + Var_Al_RelPos   + Mean_depth , data = train_data , ntree = 1500)
# Predict probabilities
probs <- predict(rfModel, newdata = test_data, type = "response")
#probs = as.numeric(probs)
# Calculate ROC
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
threshold <- 0.5
test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)
equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)
# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)
sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print(sumCorrect )


