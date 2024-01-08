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


model2 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+I(DEPTH^2)+QUAL+Var_Al_RelPos+I(Var_Al_RelPos^2)+Ref_Al_RelPos+I(Ref_Al_RelPos^2) +Coverage + Mean_depth + I(Mean_depth^2), data = dfFilt, family = "binomial")
#summary(model2)

step_model2  = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(step_model2)

#auc
set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Ref_Al_RelPos + I(Ref_Al_RelPos^2) + Coverage + I(Mean_depth^2) , data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)

# vif testing
vif_model2 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos + Ref_Al_RelPos +Coverage + Mean_depth  , data = dfFilt, family = "binomial")
print(vif(vif_model2))

vif_model3 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos +Coverage + Mean_depth  , data = dfFilt, family = "binomial")
print(vif(vif_model3))

vif_model4 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos +Coverage  , data = dfFilt, family = "binomial")
print(vif(vif_model4))

vif_model4 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+QUAL+Var_Al_RelPos +Coverage + Mean_depth  , data = dfFilt, family = "binomial")
print(vif(vif_model4))

aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS  + QUAL + Var_Al_RelPos + I(Mean_depth^2)   , data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)


finalStepModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Coverage, data = dfFilt , family = "binomial")
step_model3  = stepAIC(finalStepModel , direction = "both" , trace = T, steps = 10000)
summary(step_model3)

# try with train data
finalStepModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Coverage, data = train_data , family = "binomial")
step_model3  = stepAIC(finalStepModel , direction = "both" , trace = T, steps = 10000)
summary(step_model3)


## Ampseq

df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))


model2 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+I(DEPTH^2)+QUAL+Var_Al_RelPos+I(Var_Al_RelPos^2)+
                Ref_Al_RelPos+I(Ref_Al_RelPos^2) +Coverage + Mean_depth + I(Mean_depth^2), data = dfFilt, family = "binomial")
#summary(model2)

step_model2  = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(step_model2)
# auc
set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) 
               + Ref_Al_RelPos + I(Ref_Al_RelPos^2) + Mean_depth + I(Mean_depth^2) , data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response") 
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE) # 0.951

# vif testing
vif_model2 <- glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + QUAL + Var_Al_RelPos +
                   Ref_Al_RelPos + Mean_depth   , data = dfFilt, family = "binomial")
print(vif(vif_model2))

vif_model3 <- glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + QUAL + Var_Al_RelPos +
                     Mean_depth   , data = dfFilt, family = "binomial")
print(vif(vif_model3))

vif_model4 <- glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH  + Var_Al_RelPos + Mean_depth    , data = dfFilt, family = "binomial")
print(vif(vif_model4))


aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos  +
                 Mean_depth  , data = train_data , family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response") 
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)


finalStepModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) +
                Mean_depth + I(Mean_depth^2), data = train_data , family = "binomial")
step_model3  = stepAIC(finalStepModel , direction = "both" , trace = T, steps = 10000)
summary(step_model3)


finalStepModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos + I(Var_Al_RelPos^2) +
                       Mean_depth + I(Mean_depth^2), data = dfFilt , family = "binomial")
step_model3  = stepAIC(finalStepModel , direction = "both" , trace = T, steps = 10000)
summary(step_model3)
