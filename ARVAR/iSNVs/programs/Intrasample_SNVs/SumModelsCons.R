# make a function to test individual variables with spearman rank correlation.


oneVarCoef = function(dfPath) {
  varsList = c("ALLELE.FREQUENCY","STRAND.BIAS", "DEPTH", "QUAL", "Var_Al_RelPos", "Ref_Al_RelPos", "coverage", "meandepth", "meanbaseq", "meanmapq")
  varNames = c("Allele Frequency at SNV position", "Strand bias at SNV position", "Sequencing depth at SNV position", "Base quality at SNV position", 
               "Mean realtive position for variant base", "Mean relative position for reference base", "percentage of covered genome per sample", "mean seuqencing depth per sample",
               'mean base quality per sample', "mean mapping quality per sample")
  df = read.csv(dfPath)
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
  dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
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
  combData$Full_variable_name = varNames
  return(combData)
}

multivarModel = function(dfPath) {
  varNames = c('Intercept',"Allele Frequency at SNV position", "Strand bias at SNV position", "Base quality at SNV position", 
               "Mean realtive position for variant base", "Mean relative position for reference base", "mean seuqencing depth per sample",
               'mean base quality per sample')
  df = read.csv(dfPath)
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  dfFilt =   dfFilt[!is.na(dfFilt$Var_Al_RelPos),]
  dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
  curModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + QUAL + Var_Al_RelPos + Ref_Al_RelPos  + meandepth + meanbaseq, data = dfFilt, family = "binomial")
  curSummary= data.frame(summary(curModel)$coefficients)
  curSummary$Variable = rownames(curSummary)
  colnames(curSummary)[4] = "Pvalue"
  curSummary$Full_variable_name = varNames
  return(curSummary)
}


ampseq = oneVarCoef("IntraSnv_results/ampseq_ConsTest_freq_1_0.csv")
ampseqMulti = multivarModel("IntraSnv_results/ampseq_ConsTest_freq_1_0.csv")

ampseqFreq = oneVarCoef("IntraSnv_results/ampseq_ConsTest_freq_1_0.01.csv")
ampseqMultiFreq = multivarModel("IntraSnv_results/ampseq_ConsTest_freq_1_0.01.csv")


metaseq = oneVarCoef("IntraSnv_results/metaseq_ConsTest_freq_1_0.csv")
metaseqMulti = multivarModel("IntraSnv_results/metaseq_ConsTest_freq_1_0.csv")

metaseqFreq = oneVarCoef("IntraSnv_results/metaseq_ConsTest_freq_1_0.01.csv")
metaseqMultiFreq = multivarModel("IntraSnv_results/metaseq_ConsTest_freq_1_0.01.csv")

targDir = "IntraSnv_results/SumModels/"
dir.create(targDir, recursive = T, showWarnings = F)

write.csv(ampseq, paste0(targDir, "ampseq_univar.csv"), row.names = F)
write.csv(ampseqMulti, paste0(targDir, "ampseq_multivar.csv"), row.names = F)
write.csv(ampseqFreq, paste0(targDir, "ampseq_univar_freq_1_0.01.csv"), row.names = F)
write.csv(ampseqMultiFreq, paste0(targDir, "ampseq_multivar_freq_1_0.01.csv"), row.names = F)

write.csv(metaseq, paste0(targDir, "metaseq_univar.csv"), row.names = F)
write.csv(metaseqMulti, paste0(targDir, "metaseq_multivar.csv"), row.names = F)
write.csv(metaseqFreq, paste0(targDir, "metaseq_univar_freq_1_0.01.csv"), row.names = F)
write.csv(metaseqMultiFreq, paste0(targDir, "metaseq_multivar_freq_1_0.01.csv"), row.names = F)

system("aws s3 cp --recursive IntraSnv_results/SumModels/ s3://abombin/ARVAR/iSNVs/IntraSnv_results/SumModels/")


system("aws s3 sync IntraSnv_results/SumModels/ s3://abombin/ARVAR/iSNVs/IntraSnv_results/SumModels/")

system("aws s3 cp --recursive IntraSnv_results/class_acuracy/ s3://abombin/ARVAR/iSNVs/IntraSnv_results/class_acuracy/")