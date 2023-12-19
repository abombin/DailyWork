library(randomForest)

set.seed(42)

dfFilt = df[!is.na(df$Var_Al_RelPos),]

set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

# Fit a Random Forest model
#rfModel <- randomForest(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos+Ref_Al_RelPos+ +Coverage + Mean_depth , data = train_data , ntree = 1500)  # You can adjust the number of trees (ntree) as needed
rfModel <- randomForest(ConsTest  ~ ALLELE.FREQUENCY + STRAND.BIAS  + QUAL + Var_Al_RelPos   + Mean_depth , data = train_data , ntree = 1500)


rfModel <- randomForest(ConsTest  ~ ALLELE.FREQUENCY   + QUAL + Var_Al_RelPos   + Mean_depth , data = train_data , ntree = 1500)

# Predict probabilities
probs <- predict(rfModel, newdata = test_data, type = "response")

# Calculate ROC
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)

threshold <- 0.5
test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)
equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)
# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)
sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print(sumCorrect )