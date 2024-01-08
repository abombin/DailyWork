library(MASS)
library(pROC)
library(caret)

getAverageAuc = function(DF, reps, protocol) {
  df = read.csv(DF)
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
  dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
  combAuc = numeric()
  for ( i in 1:reps ) {
    curSeed = (i + 8) * 4
    print(curSeed)
    set.seed(curSeed)
    train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
    train_data <- dfFilt[train_idx, ]
    test_data <- dfFilt[-train_idx, ]
    
    if (protocol == "metaseq") {
      aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + I(DEPTH^2) + QUAL + Var_Al_RelPos + I(Var_Al_RelPos^2) + Ref_Al_RelPos + I(Ref_Al_RelPos^2), data = train_data , family = "binomial")
    }
    else if (protocol == "ampseq") {
      aucModel = glm(ConsTest ~ ALLELE.FREQUENCY + DEPTH + Var_Al_RelPos + I(Var_Al_RelPos^2) + Mean_depth + I(Mean_depth^2), data = train_data , family = "binomial")
    }
    probs <- predict(aucModel, newdata = test_data, type = "response")
    roc_obj <- roc(test_data$ConsTest ~ probs, plot = F, print.auc = F)
    curAuc = as.numeric(roc_obj$auc)
    combAuc= c(combAuc, curAuc)
  }
  meanAuc = mean(combAuc)
  return(meanAuc)
}

metaseqAuc = getAverageAuc(DF = "snvs_comb_res/metaseq_ampseq_overlap_comb_derep_decont_covFilt_97.csv", reps = 100, protocol = "metaseq")

ampseqAuc = getAverageAuc(DF = "snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv", reps = 100, protocol = "ampseq")