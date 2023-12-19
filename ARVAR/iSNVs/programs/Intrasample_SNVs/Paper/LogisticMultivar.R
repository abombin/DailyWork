dir.create("Paper/Tables", recursive = T, showWarnings = F)

df = read.csv("IntraSnv_results/ampseq_ConsTest_freq_1_0.csv")
#df = read.csv("IntraSnv_results/metaseq_ConsTest_freq_1_0.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))


finalModel = glm(ConsTest ~ALLELE.FREQUENCY + STRAND.BIAS + QUAL +Var_Al_RelPos + Ref_Al_RelPos  + meandepth + meanbaseq, data = dfFilt, family = "binomial")
summary(finalModel )


varNames = c('Intercept',"Allele Frequency at SNV position", "Strand bias at SNV position", "Base quality at SNV position", 
             "Mean realtive position for variant base", "Mean relative position for reference base", "mean seuqencing depth per sample",
             'mean base quality per sample')
curSummary= data.frame(summary(finalModel)$coefficients)
curSummary$Variable = rownames(curSummary)
colnames(curSummary)[4] = "Pvalue"
curSummary$Full_variable_name = varNames

curSummary[order(abs(curSummary$Estimate), decreasing = T),]


#write.csv(curSummary, "Paper/Tables/Ampseq_Multivar_LogisticCoef_Table.csv", row.names = F)

#write.csv(curSummary, "Paper/Tables/Metaseq_Multivar_LogisticCoef_Table.csv", row.names = F)