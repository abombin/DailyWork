# get AUC
library(pROC)
library(caret)

df = read.csv("snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv")

#df = read.csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv")

df1 <- df[complete.cases(df[, c("Var_Al_RelPos", "Ref_Al_RelPos")]), ]

set.seed(42)
train_idx <- createDataPartition(df$ConsTest, p = 0.7, list = FALSE)
train_data <- df[train_idx, ]
test_data <- df[-train_idx, ]



multivarModel =  glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos + Ref_Al_RelPos, data = train_data , family = "binomial")
summary(multivarModel)

probs <- predict(multivarModel, newdata = test_data, type = "response")

roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)