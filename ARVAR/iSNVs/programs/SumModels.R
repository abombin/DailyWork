library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)


## Metaseq
df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

model1 <- glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Coverage, data = dfFilt, family = "binomial")
summary(model1)

step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model)

sumModel = function(model) {
  modSum = summary(model)
  modRes = data.frame(modSum[["coefficients"]])
  colnames(modRes)[4] = "Pval"
  modRes$Variable = rownames(modRes)
  modRes$AIC = modSum$aic
  return(modRes)
}

multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Coverage ,data = dfFilt, family = "binomial")

modRes2 = sumModel(multivarModel)

# test variance inflation score 
print(vif(multivarModel))
vif_model <- glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos  + Coverage, data = dfFilt, family = "binomial")
print(vif(vif_model))
cur_vif = data.frame(vif(vif_model))
colnames(cur_vif)[1] = "VIF"

#AUC
set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]
aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Ref_Al_RelPos + I(Ref_Al_RelPos^2), data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)

#
dir.create("model_summary", showWarnings = F)
write.csv(modRes2, "model_summary/metaseq_ampseq_97_bestAic_sqBiomodel_sum.csv", row.names = F)
write.csv(cur_vif, "model_summary/metaseq_ampseq_97_bestAic_sqBiomodel_VIF.csv", row.names = F)

png(file = "model_summary/metaseq_ampseq_97_bestAic_sqBiomodel_auc.png" ,height = 10, width = 10, units = "in", res = 300)
print(roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE))
dev.off()


## Ampseq
df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

finalStepModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos + I(Var_Al_RelPos^2) +
                       Mean_depth + I(Mean_depth^2), data = dfFilt , family = "binomial")
step_model3  = stepAIC(finalStepModel , direction = "both" , trace = T, steps = 10000)
summary(step_model3)

modRes2 = sumModel(finalStepModel)


# test variance inflation score 
print(vif(finalStepModel))
vif_model <- glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos +
                   Mean_depth , data = dfFilt, family = "binomial")
print(vif(vif_model))
cur_vif = data.frame(vif(vif_model))
colnames(cur_vif)[1] = "VIF"

#AUC
set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]
aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos + I(Var_Al_RelPos^2) +
                 Mean_depth + I(Mean_depth^2), data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)


write.csv(modRes2, "model_summary/ampseq_metaseq_97_bestAic_sqBiomodel_sum.csv", row.names = F)
write.csv(cur_vif, "model_summary/ampseq_metaseq_97_bestAic_sqBiomodel_VIF.csv", row.names = F)

png(file = "model_summary/ampseq_metaseq_97_bestAic_sqBiomodel_auc.png" ,height = 10, width = 10, units = "in", res = 300)
print(roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE))
dev.off()