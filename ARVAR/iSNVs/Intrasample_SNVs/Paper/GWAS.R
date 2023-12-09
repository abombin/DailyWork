# gwas
library(logistf)

df = read.csv("Paper/Additional_analysis/GWAS/SnvsMatMaf0.01_Wu.csv")

cols= colnames(df)
samples = cols[4:length(cols)]
samples = data.frame(gsub("_", "-", samples))
colnames(samples) = "Sample"

# metadata
metadat = read.csv("Final_vaxbt_dataset_AP_metadata.csv")
metadatFilt = unique(metadat[, c("AP_lab_id", "Ct_value", "WHO_variant", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity" )])
metadatFilt[metadatFilt == ""] <- NA
metadatFilt[metadatFilt == "."] <- NA
metadatFilt$Sample1 = gsub('_', "-", metadatFilt$AP_lab_id)
metadatFilt$Sample = gsub('_', "-", metadatFilt$AP_lab_id)

combMeta = plyr::join(samples, metadatFilt, by = "Sample", type = "left", match = 'all')
identical(combMeta$Sample, samples$Sample)

combMeta$Vaccinated = NA
combMeta$Vaccinated[combMeta$vax_doses_received == 0]<-"No"
combMeta$Vaccinated[combMeta$vax_doses_received > 0]<-"Yes"

i = 1

# firth regression
combDat = data.frame()
for ( i in 1:nrow(df) ) {
  curSNV = paste(df$POSITION[i], df$REF.NT[i], df$VAR.NT[i], sep = "__")
  curVals = data.frame(t(df[i, 4:ncol(df)]))
  colnames(curVals) = "Values"
  curVals$Sample = rownames(curVals)
  curVals$Sample = gsub("_", "-", curVals$Sample)
  testDf = plyr::join(curVals, combMeta, by = "Sample", type = "left", match = 'all')
  lf <- logistf(formula = Values ~ Vaccinated, data =  testDf)
  sumLf = summary(lf)
  coefLf = data.frame(cbind(sumLf$coefficients, sumLf$prob))
  lm = glm(Values ~ Vaccinated, data =  testDf, family = binomial)
  summary(lm)
  colnames(coefLf) = c("Coef", "Pval")
  coefLf$Pos_Ref_Var = curSNV
  coefLf$Vars = rownames(coefLf)
  coefLf = coefLf[coefLf$Vars != "(Intercept)",]
  combDat = rbind(combDat, coefLf)
}
combDat$bonf_p = p.adjust(combDat$Pval, method = "bonferroni")
combDat$fdr_p = p.adjust(combDat$Pval, method = "fdr")

# logistic regression
combDat = data.frame()
for ( i in 1:nrow(df) ) {
  curSNV = paste(df$POSITION[i], df$REF.NT[i], df$VAR.NT[i], sep = "__")
  curVals = data.frame(t(df[i, 4:ncol(df)]))
  colnames(curVals) = "Values"
  curVals$Sample = rownames(curVals)
  curVals$Sample = gsub("_", "-", curVals$Sample)
  testDf = plyr::join(curVals, combMeta, by = "Sample", type = "left", match = 'all')
  lf = glm(Values ~ Vaccinated, data =  testDf, family = binomial)
  sumLf = summary(lf)
  coefLf = data.frame(sumLf$coefficients)
  colnames(coefLf)[4] = "Pval"
  coefLf$Pos_Ref_Var = curSNV
  coefLf$Vars = rownames(coefLf)
  coefLf = coefLf[coefLf$Vars != "(Intercept)",]
  combDat = rbind(combDat, coefLf)
}
combDat$bonf_p = p.adjust(combDat$Pval, method = "bonferroni")
combDat$fdr_p = p.adjust(combDat$Pval, method = "fdr")