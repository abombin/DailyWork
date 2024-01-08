# for this code instead of using treshold of 0.5 need to find the best treshold based on F1 score


DF = "snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv"

curSeed = 42

df = read.csv(DF)
dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

set.seed(curSeed)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Ref_Al_RelPos + I(Ref_Al_RelPos^2), data = train_data , family = "binomial")

probs <- predict(aucModel, newdata = test_data, type = "response")

threshold <- 0.5
test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)

equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)

# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)

sumCorrect = equal_count / (equal_count + not_equal_count) * 100

roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)

## ampseq
DF = "snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv"

curSeed = 42

df = read.csv(DF)
dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

set.seed(curSeed)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos + I(Var_Al_RelPos^2) + Mean_depth + I(Mean_depth^2), data = train_data , family = "binomial")


probs <- predict(aucModel, newdata = test_data, type = "response")

threshold <- 0.5
test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)

equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)

# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)

sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print(sumCorrect)

roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)