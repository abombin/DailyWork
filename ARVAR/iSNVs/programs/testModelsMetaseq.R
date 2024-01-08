library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)


sumModel = function(model) {
  modSum = summary(model)
  modRes = data.frame(modSum[["coefficients"]])
  colnames(modRes)[4] = "Pval"
  modRes$Variable = rownames(modRes)
  modRes$AIC = modSum$aic
  return(modRes)
}


## Metaseq
df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt = dfFilt[dfFilt$ALLELE.FREQUENCY <=0.98,]

dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))


model2 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+I(DEPTH^2)+QUAL+Var_Al_RelPos+I(Var_Al_RelPos^2)+Ref_Al_RelPos+I(Ref_Al_RelPos^2) +Coverage + Mean_depth + I(Mean_depth^2), data = dfFilt, family = "binomial")
#summary(model2)

step_model2  = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(step_model2)

#auc
set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + 
                 QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Ref_Al_RelPos + 
                 I(Ref_Al_RelPos^2) + Coverage + Mean_depth + I(Mean_depth^2), data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)


# vif testing
vif_model2 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos + Ref_Al_RelPos +Coverage + Mean_depth  , data = dfFilt, family = "binomial")
print(vif(vif_model2))

# vif testing
vif_model2 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+QUAL+Var_Al_RelPos + Ref_Al_RelPos +Coverage + Mean_depth  , data = dfFilt, family = "binomial")
print(vif(vif_model2))


set.seed(42)
aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS  + 
                 QUAL + Var_Al_RelPos   + I(Mean_depth^2) , data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)


#####

threshold <- 0.5
test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)

equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)

# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)

sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print(sumCorrect )

# filter test data
testDatFilt = test_data[test_data$ALLELE.FREQUENCY <= 0.98, ]
testDatFilt = test_data[test_data$Freq_adj <= 0.98, ]
equal_count <- sum(testDatFilt$ConsTest == testDatFilt$PredictedConsTest)
# Count the number of times values in Column1 are NOT equal to Column2
not_equal_count <- sum(testDatFilt$ConsTest != testDatFilt$PredictedConsTest)
sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print(sumCorrect )

## save results

finalModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS  + 
                 QUAL + Var_Al_RelPos  + I(Var_Al_RelPos^2) +  I(Mean_depth^2)  , data = dfFilt , family = "binomial")
modelSum = sumModel(finalModel)

final_vif <- glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS  + 
                   QUAL + Var_Al_RelPos + I(Mean_depth^2)  , data = dfFilt, family = "binomial")
print(vif(final_vif))


write.csv(modelSum , "model_summary/metaseq_ampseq_97_freq_0.02-0.98_sum.csv", row.names = F)
#write.csv(final_vif, "model_summary/metaseq_ampseq_97_freq_0.02-0.98_VIF.csv", row.names = F)

png(file = "model_summary/metaseq_ampseq_97_freq_0.02-0.98_lin_auc.png" ,height = 10, width = 10, units = "in", res = 300)
print(roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE))
dev.off()