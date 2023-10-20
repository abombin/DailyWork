library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)


df = read.csv("IntraSnv_results/metaseq_ConsTest_freq_1_0.01.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

#model1 <- glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + I(DEPTH^2)  + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Ref_Al_RelPos + coverage + meandepth + meanbaseq +meanmapq, data = dfFilt, family = "binomial")
model1 <- glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos + Ref_Al_RelPos + coverage + meandepth + meanbaseq +meanmapq, data = dfFilt, family = "binomial")
summary(model1)

oneVarCoef = function(varsList, dfFilt) {
  combData = data.frame()
  for (var in varsList) {
    curForm = as.formula(paste0("ConsTest ~ ", var))
    curModel = glm(curForm, data = dfFilt, family = "binomial")
    curSummary = data.frame(summary(curModel)$coefficients)
    curSummary = curSummary[2,]
    curSummary$Variable = rownames(curSummary)
    combData= rbind(combData, curSummary)
    
  }
  colnames(combData)[4] = "Pvalue"
  return(combData)
}

varsList = c("ALLELE.FREQUENCY","STRAND.BIAS", "DEPTH", "QUAL", "Var_Al_RelPos", "Ref_Al_RelPos", "coverage", "meandepth", "meanbaseq", "meanmapq")
sumVars = oneVarCoef(varsList, dfFilt)

step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model)

# best linear model ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos + Ref_Al_RelPos + coverage + meandepth + meanbaseq + meanmapq

vif_model <- glm(ConsTest ~  ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos + Ref_Al_RelPos + coverage + meandepth + meanbaseq + meanmapq, data = dfFilt, family = "binomial")
print(vif(vif_model))

vif_model <- glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS  + QUAL + Var_Al_RelPos + Ref_Al_RelPos + coverage + meandepth + meanbaseq, data = dfFilt, family = "binomial")
print(vif(vif_model))

### AUC

set.seed(42)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]

aucModel <- glm(ConsTest ~   ALLELE.FREQUENCY + STRAND.BIAS  + QUAL + Var_Al_RelPos + Ref_Al_RelPos + coverage + meandepth + meanbaseq, data = dfFilt, family = "binomial")
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
print(roc_obj$auc)

aucModel <- randomForest(formula = ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + QUAL + Var_Al_RelPos + Ref_Al_RelPos  + meandepth + meanbaseq, data = train_data, ntree = 500)
probs <- predict(aucModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
print(roc_obj$auc)
importance_scores <- importance(aucModel)
importance_scores


finalModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + QUAL + Var_Al_RelPos + Ref_Al_RelPos  + meandepth + meanbaseq, data = dfFilt, family = "binomial")
summary(finalModel )

