library(rstatix)
library(dplyr)
library(car)
library(MASS)

minFreq = 0.01
addMiss = T


ShPi = function(targDf, freqCol) {
  combPi = numeric()
  for ( i in 1:nrow(targDf) ) {
    curFreq = targDf[i, ][[freqCol]]
    if ( curFreq == 0) {
      pi = 0
    } else {
      pi = curFreq * log(curFreq)
    }
    combPi = c(combPi, pi)
  }
  return(combPi)
}

getShannon = function(df) {
  allSamples = character()
  allShannon = numeric()
  allDepth = numeric()
  samplesList = unique(df$Sample)
  for (curSample in samplesList) {
    curDf = df[df$Sample == curSample,]
    curDepth = unique(curDf$meandepth)
    curShannon = -1 * sum(curDf$Pi.Ln.Pi.)
    allSamples = c(allSamples, curSample)
    allShannon = c(allShannon, curShannon)
    allDepth = c(allDepth, curDepth)
  }
  combDat = data.frame(allSamples, allShannon, allDepth)
  colnames(combDat) = c("Sample", "Shannon", "meandepth")
  return(combDat)
}


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

# metadata

metadat = read.csv("Final_vaxbt_dataset_AP_metadata.csv")
metadatFilt = unique(metadat[, c("AP_lab_id", "Ct_value", "WHO_variant", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity" )])
metadatFilt[metadatFilt == ""] <- NA
metadatFilt[metadatFilt == "."] <- NA
metadatFilt$Sample1 = gsub('_', "-", metadatFilt$AP_lab_id)
metadatFilt$Sample = gsub('_', "-", metadatFilt$AP_lab_id)

# get shannon
metaseq = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions.csv')
ampseq = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv')

# filter by frequency
metaseqOverlap = metaseq[metaseq$ALLELE.FREQUENCY >= minFreq, ]
ampseqOverlap = ampseq[ampseq$ALLELE.FREQUENCY >= minFreq, ]

metaseqOverlap$Pi.Ln.Pi. = ShPi(targDf=metaseqOverlap, freqCol="ALLELE.FREQUENCY")
ampseqOverlap$Pi.Ln.Pi. = ShPi(targDf=ampseqOverlap, freqCol="ALLELE.FREQUENCY")

# calculate shannon diversity 
metaseqOverlapShan = getShannon(df=metaseqOverlap)
ampseqOverlapShan = getShannon(df = ampseqOverlap)

# add samples that were filtered out during SNVs predictions

ampseqAllSamples = read.csv('IntraSnv_results/ampseq_comb_derep.csv')
ampseqAllSamples = ampseqAllSamples[ampseqAllSamples$coverage >= 97,]
ampseqAllSamples = unique(ampseqAllSamples[, c("Sample",  "meandepth")])
metaseqAllSamples = read.csv('IntraSnv_results/metaseq_comb_derep.csv')
metaseqAllSamples = metaseqAllSamples[metaseqAllSamples$coverage >= 97,]
metaseqAllSamples = unique(metaseqAllSamples[, c("Sample",  "meandepth")])

#min(ampseqAllSamples$coverage)
#min(metaseqAllSamples$coverage)


addMissing = function(df, refDf) {
  combSamples = character()
  combDepth = numeric()
  combShannon = numeric()
  samplesList = unique(df$Sample)
  for ( i in 1:nrow(refDf) ) {
    curSample = refDf$Sample[i]
    if (!curSample%in%samplesList) {
      curDepth = refDf$meandepth[i]
      curShannon = 0 
      
      combSamples = c(combSamples, curSample)
      combDepth = c(combDepth, curDepth)
      combShannon = c(combShannon, curShannon)
    }
  }
  missingDat = data.frame(Sample = combSamples, Shannon = combShannon, meandepth = combDepth)
  print(nrow(missingDat))
  combDf = rbind(df, missingDat)
  return(combDf)
}

# ADD missing
if (addMiss == T) {
  metaseqOverlapShan = addMissing(df=metaseqOverlapShan, refDf=metaseqAllSamples)
  ampseqOverlapShan = addMissing(df=ampseqOverlapShan, refDf=ampseqAllSamples)
}


# add metadata
metaCombShan = plyr::join(metaseqOverlapShan, metadatFilt, by = "Sample", type = "left", match = 'all')
ampCombShan = plyr::join(ampseqOverlapShan, metadatFilt, by = "Sample", type = "left", match = 'all')

#metaCombShan = metaCombShan[!is.na(metaCombShan$AP_lab_id),]
#ampCombShan = ampCombShan[!is.na(ampCombShan$AP_lab_id),]

columns_to_check <- c("Ct_value", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
metaCombShan[columns_to_check] <- lapply(metaCombShan[columns_to_check], as.numeric)
ampCombShan[columns_to_check] <- lapply(ampCombShan[columns_to_check], as.numeric)

metaCombShan$Ct_depth_adj = metaCombShan$Ct_value / metaCombShan$meandepth
ampCombShan$Ct_depth_adj = ampCombShan$Ct_value / ampCombShan$meandepth

# add vaccination status
metaCombShan$Vaccinated = NA
metaCombShan$Vaccinated[metaCombShan$vax_doses_received == 0]<-"No"
metaCombShan$Vaccinated[metaCombShan$vax_doses_received > 0]<-"Yes"

ampCombShan$Vaccinated = NA
ampCombShan$Vaccinated[ampCombShan$vax_doses_received == 0]<-"No"
ampCombShan$Vaccinated[ampCombShan$vax_doses_received > 0]<-"Yes"


# evaluate all samples together
metaCombShan$Protocol = "metaseq"
ampCombShan$Protocol = "ampseq"
combShannon = rbind(metaCombShan, ampCombShan)

##
combShannon$NormShan = combShannon$Shannon
minVal = min(combShannon$Shannon [combShannon$Shannon != 0]) * 0.5
combShannon$NormShan[combShannon$NormShan == 0] = minVal
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))

# univariate nin parametric tests 

combShannon = combShannon[!combShannon$days_post_symptom_onset > 100,]
varList = c("Ct_depth_adj", "Ct_value","vax_doses_received", "meandepth", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
varList = c("vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
ShanCorAll = runSpearman(df=combShannon, varList=varList, testVar="Shannon", test = "spearman")
wilcVax_2_All = runWilcox(curPred ="Vaccinated", combShannon= combShannon)
wilcVax = runWilcox(curPred ="vax_doses_received", combShannon= combShannon)
combShannon = combShannon[combShannon$WHO_variant != "Beta",]
wilcVarAll = runWilcox(curPred ="WHO_variant", combShannon= combShannon)


# multivariate model
shapiro.test(combShannon$NormShan)
ks.test(combShannon$NormShan, "pnorm")

DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
varSel = varSel[!varSel$days_post_symptom_onset > 100,]
varSel = drop_na(varSel)
multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
model2 = lm(NormShan~Protocol + meandepth +WHO_variant+Vaccinated + days_post_symptom_onset , data = varSel) 

# residuals <- residuals(model2)
# ks.test(residuals, "pnorm")
# shapiro.test(residuals)

fullMult = summary(model2)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns
fullMult$bonferroni_p = p.adjust(fullMult$Pval, method = "bonferroni")
fullMult$fdr_p = p.adjust(fullMult$Pval, method = "fdr")
fullMult$Predictor = rownames(fullMult)
predictorNames = c("Intercept","Metagenomic Sequencing","Omicron Variant ","Vaccinated","Days post symptom onset","Days post symptom onset^2")
fullMult$Predictor = predictorNames
#write.csv(fullMult, "Paper/Tables/Original_Multivar.csv", row.names = F)


# step AIC
model1 = lm(NormShan~Protocol+meandepth + Ct_depth_adj +WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) , data = varSel) 
#model1 = lm(NormShan~Protocol+meandepth+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) , data = varSel) 
#model1 = lm(NormShan~Protocol+meandepth + Ct_value+WHO_variant+Vaccinated + days_post_symptom_onset  , data = varSel) 
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model)
residuals <- residuals(step_model)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
fullMult = summary(step_model)
fullMult = data.frame(fullMult$coefficients)
multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
colnames(fullMult) = multColumns
# write.csv(fullMult, "Paper/Tables/shannon_Multivar_Freq_1-0.02_addMiss.csv", row.names = T)

fullMult

if (addMiss == T) {
  write.csv(ShanCorAll, paste0("Paper/Tables/Shannon_Spearman_Freq_1-",minFreq, "_addMiss.csv"), row.names = F)
  write.csv(wilcVax_2_All, paste0("Paper/Tables/Shannon_Wilcox_Vax_Freq_1-",minFreq, "_addMiss.csv"), row.names = F)
  write.csv(wilcVax, paste0("Paper/Tables/Shannon_Wilcox_VaxGroups_Freq_1-", minFreq, "_addMiss.csv"), row.names = F)
  write.csv(wilcVarAll, paste0("Paper/Tables/Shannon_Wilcox_WHO_Freq_1-", minFreq, "_addMiss.csv"),  row.names = F)
  
  write.csv(fullMult, paste0("Paper/Tables/shannon_Multivar_Freq_1-", minFreq, "_addMiss.csv"), row.names = T)
} else {
  write.csv(ShanCorAll, paste0("Paper/Tables/Shannon_Spearman_Freq_1-",minFreq, ".csv"), row.names = F)
  write.csv(wilcVax_2_All, paste0("Paper/Tables/Shannon_Wilcox_Vax_Freq_1-",minFreq, ".csv"), row.names = F)
  write.csv(wilcVax, paste0("Paper/Tables/Shannon_Wilcox_VaxGroups_Freq_1-", minFreq, ".csv"), row.names = F)
  write.csv(wilcVarAll, paste0("Paper/Tables/Shannon_Wilcox_WHO_Freq_1-", minFreq, ".csv"),  row.names = F)
  
  write.csv(fullMult, paste0("Paper/Tables/shannon_Multivar_Freq_1-", minFreq, ".csv"), row.names = T)
}



system("aws s3 sync Paper/Tables/ s3://abombin/ARVAR/iSNVs/Paper/Tables/")

shapiro.test(varSel$NormShan)
ks.test(varSel$NormShan, "pnorm")

# deleted Original_Multivar as accidentaly substututed it with the wrong one