combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")

runSpearman = function(df, varList, testVar, test) {
  Variable = character()
  Estimate = numeric()
  Pval = numeric()
  curResp = df[[testVar]]
  for (i in varList ) {
    #print(i)
    curPredict = df[[i]]
    curCor = cor.test(curResp, curPredict, method = test)
    cur_cor_coef = curCor$estimate
    cur_pval = curCor$p.value
    Variable = c(Variable, i)
    Estimate = c(Estimate, cur_cor_coef)
    Pval = c(Pval, cur_pval)
  }
  combDat = data.frame(Variable, Estimate, Pval)
}

runWilcox = function(curPred, combShannon) {
  wilcForm = as.formula(paste("Shannon ~", curPred ))
  wilcox_res<-combShannon %>% rstatix::pairwise_wilcox_test(wilcForm, p.adjust.method = "fdr")
  sumRes = result <- aggregate(wilcForm, data = combShannon, FUN = mean)
  colnames(sumRes)[1] = "Group"
  combDiff = numeric()
  for (i in 1:nrow(wilcox_res)) {
    group1 =  wilcox_res$group1[i]
    group2 = wilcox_res$group2[i]
    curDiff = sumRes$Shannon[sumRes$Group == group1] - sumRes$Shannon[sumRes$Group == group2]
    combDiff= c(combDiff, curDiff)
  }
  wilcox_res$Groups_difference = combDiff
  wilcox_res = wilcox_res[wilcox_res$n1 >= 4 & wilcox_res$n2 >= 4, ]
  colnames(wilcox_res)[1] = "Variable"
  
  sel_columns = c("Variable", "group1", "group2" , "p", "Groups_difference")
  wilcox_res =  wilcox_res[, sel_columns]
  
  return(wilcox_res)
}

#varList = c("Ct_depth_adj", "Ct_value","vax_doses_received", "meandepth", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
varList = c("vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
ShanCorAll = runSpearman(df=combShannon, varList=varList, testVar="Shannon", test = "spearman")
wilcVax_2_All = runWilcox(curPred ="Vaccinated", combShannon= combShannon)
wilcVarAll = runWilcox(curPred ="WHO_variant", combShannon= combShannon)

write.csv(ShanCorAll, "Paper/Tables/Shannon_Spearman.csv", row.names = F)
write.csv(wilcVax_2_All, "Paper/Tables/Shannon_Wilcox_Vax.csv", row.names = F)
write.csv(wilcVarAll, "Paper/Tables/Shannon_Wilcox_WHO.csv",  row.names = F)
