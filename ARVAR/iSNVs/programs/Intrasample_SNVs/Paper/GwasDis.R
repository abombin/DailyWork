
library(ordinal)

setwd("/home/ubuntu/extraVol/ARVAR/iSNVs/Paper/Additional_analysis/GWAS")

df = read.csv("SnvsMatMaf0.01_MinFreq0.01_Wu_oversamp.csv")

protDat = read.csv("AmpMetaCombSnvsPostPredDedup_MinFreq0.01_Wu_oversamp.csv")
protDat = unique(protDat[, c("Sample", "Protocol")])

cols= colnames(df)
samples = cols[4:length(cols)]
samples = data.frame(gsub("_", "-", samples))
colnames(samples) = "Sample"

metadat = read.csv("../../../Final_vaxbt_dataset_AP_metadata.csv")
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

mdsComp = read.csv("Mds_covar_Wu.csv")

combMeta = plyr::join(combMeta, mdsComp, by = "Sample", type = "left", match = 'all')
finalMeta = plyr::join(combMeta, protDat, by = "Sample", type = "left", match = 'all')
nrow(combMeta) == nrow(finalMeta)

#finalMeta  = finalMeta[finalMeta$vax_doses_received < 3,]

runGwas = function(df, combMeta, adjPop) {
  # firth regression
  combDat = data.frame()
  for ( i in 1:nrow(df) ) {
    curSNV = paste(df$POSITION[i], df$REF.NT[i], df$VAR.NT[i], sep = "__")
    curVals = data.frame(t(df[i, 4:ncol(df)]))
    colnames(curVals) = "Values"
    curVals$Sample = rownames(curVals)
    curVals$Sample = gsub("_", "-", curVals$Sample)
    testDf = plyr::join(curVals, combMeta, by = "Sample", type = "left", match = 'all')
    testDf = testDf [!is.na(testDf$disease_severity),]
    testDf$disease_severity = as.factor(as.character(testDf$disease_severity))
    if (adjPop == T) {
      lf <- clm(formula = disease_severity ~ Protocol+ X1+X2+X3+X4+X5 + Values, data =  testDf)
    } else {
      lf <- clm(formula = disease_severity ~ Protocol + Values, data =  testDf)
    }
    sumLf = summary(lf)
    coefLf = data.frame(sumLf$coefficients)
    colnames(coefLf)[4] = c("Pval")
    coefLf$Pos_Ref_Var = curSNV
    coefLf$Vars = rownames(coefLf)
    coefLf = coefLf[coefLf$Vars == "Values",]
    combDat = rbind(combDat, coefLf)
  }
  combDat$bonf_p = p.adjust(combDat$Pval, method = "bonferroni")
  combDat$fdr_p = p.adjust(combDat$Pval, method = "fdr")
  return(combDat)
}



SnvVax = runGwas(df=df, combMeta=finalMeta, adjPop=F)
SnvVaxPopAdj = runGwas(df=df, combMeta=finalMeta, adjPop=T)
SnvVaxSign = SnvVax[SnvVax$Pval < 0.05,]
SnvVaxPopAdjSign = SnvVaxPopAdj[SnvVaxPopAdj$Pval < 0.05, ]

# try with logistic regression
# SnvVax_1 = runGwasLog(df=df, combMeta=combMeta, adjPop=F)
# SnvVaxPopAdj_1 = runGwasLog(df=df, combMeta=combMeta, adjPop=T)
# SnvVaxSign_1 = SnvVax_1[SnvVax_1$Pval < 0.05,]
# SnvVaxPopAdjSign_1 = SnvVaxPopAdj_1[SnvVaxPopAdj_1$Pval < 0.05, ]

#dir.create("Results/")
write.csv(SnvVaxSign, "Results/Snv_disease_severity_maf0.01_MinFreq0.01_Wu_2023-12-28_oversamp.csv", row.names = F)
write.csv(SnvVaxPopAdjSign, "Results/Snv_disease_severity_PopAdjust_maf0.01_MinFreq0.01_Wu_2023-12-28_oversamp.csv", row.names = F)

system("aws s3 sync Results/ s3://abombin/ARVAR/iSNVs/Paper/Tables/GWAS/")


