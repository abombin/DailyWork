library(rstatix)
library(dplyr)
library(car)
library(MASS)

setwd("/home/ubuntu/extraVol/ARVAR/iSNVs/")

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
  allCoverage = numeric()
  allOrigin = character()
  allEnd = numeric()
  samplesList = unique(df$Sample)
  for (curSample in samplesList) {
    curDf = df[df$Sample == curSample,]
    curDepth = unique(curDf$meandepth)
    curShannon = -1 * sum(curDf$Pi.Ln.Pi.)
    curCoverage = unique(curDf$coverage)
    curOrigin = unique(curDf$Origin)
    curEnd = unique(curDf$endpos)
    allSamples = c(allSamples, curSample)
    allShannon = c(allShannon, curShannon)
    allDepth = c(allDepth, curDepth)
    allCoverage = c(allCoverage, curCoverage)
    allOrigin = c(allOrigin, curOrigin)
    allEnd = c(allEnd,curEnd)
    
  }
  combDat = data.frame(allSamples, allShannon, allDepth, allCoverage, allOrigin, allEnd)
  colnames(combDat) = c("Sample", "Shannon", "meandepth", "coverage", "Origin", "endpos")
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

addMissing = function(df, refDf, overlapSamp) {
  combSamples = character()
  combDepth = numeric()
  combShannon = numeric()
  combCoverage = numeric()
  combOrigin = character()
  combEnd = numeric()
  samplesList = unique(df$Sample)
  for ( i in 1:nrow(refDf) ) {
    curSample = refDf$Sample[i]
    if (!curSample%in%samplesList) {
      curDepth = refDf$meandepth[i]
      curCoverage = refDf$coverage[i]
      curEnd = refDf$endpos[i]
      if (curSample %in% overlapSamp) {
        curOrigin = "Overlap"
      } else {
        curOrigin = "Filtered"
      }
      curShannon = 0 
      
      combSamples = c(combSamples, curSample)
      combDepth = c(combDepth, curDepth)
      combShannon = c(combShannon, curShannon)
      combCoverage = c(combCoverage, curCoverage)
      combOrigin = c(combOrigin, curOrigin)
      combEnd = c(combEnd, curEnd)
    }
  }
  missingDat = data.frame(Sample = combSamples, Shannon = combShannon, meandepth = combDepth, coverage=combCoverage, Origin = combOrigin, endpos=combEnd)
  print(nrow(missingDat))
  combDf = rbind(df, missingDat)
  return(combDf)
}

makeAllTables = function(minFreq, addMiss, writeRes) {
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
  ampseqAllSamples = unique(ampseqAllSamples[, c("Sample",  "meandepth", "coverage", "endpos")])
  metaseqAllSamples = read.csv('IntraSnv_results/metaseq_comb_derep.csv')
  metaseqAllSamples = metaseqAllSamples[metaseqAllSamples$coverage >= 97,]
  metaseqAllSamples = unique(metaseqAllSamples[, c("Sample",  "meandepth", "coverage", "endpos")])
  
  overlapSamp = unique(metaseqAllSamples$Sample[metaseqAllSamples$Sample%in%ampseqAllSamples$Sample])
  
  #min(ampseqAllSamples$coverage)
  #min(metaseqAllSamples$coverage)
  
  # ADD missing
  if (addMiss == T) {
    metaseqOverlapShan = addMissing(df=metaseqOverlapShan, refDf=metaseqAllSamples, overlapSamp=overlapSamp)
    ampseqOverlapShan = addMissing(df=ampseqOverlapShan, refDf=ampseqAllSamples, overlapSamp=overlapSamp)
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
  metaCombShan$Protocol = "Metaseq"
  ampCombShan$Protocol = "Ampseq"
  combShannon = rbind(metaCombShan, ampCombShan)
  
  
  ### Stopped here need to separate overlapping samples and predicted samples count averages and combine. 
  # create 3rd category for Protocol named Overlap
  
  
  overlapShann = combShannon[combShannon$Origin == "Overlap",]
  predictShan = combShannon[combShannon$Origin != "Overlap",]
  
  MeanShann <- aggregate(cbind(Shannon, meandepth, coverage, Ct_value, Ct_depth_adj, endpos) ~ Sample, data = overlapShann , FUN = mean)
  MeanShann$Protocol = "Overlap"
  
  colnames(combShannon)
  
  ref = unique(combShannon[, c("Sample", "WHO_variant", "vax_doses_received", "days_post_symptom_onset", "Vaccinated", "Origin", "AP_lab_id", "disease_severity", "Sample1", "days_since_last_vax" )])
  
  JoinShannon = plyr::join(MeanShann, ref, by = "Sample", type = "left", match = "all")
  nrow(MeanShann) == nrow(JoinShannon)
  
  combShannon=rbind(JoinShannon, predictShan)
  
  ##
  combShannon$NormShan = combShannon$Shannon / combShannon$endpos
  minVal = min(combShannon$Shannon [combShannon$Shannon != 0]) * 0.5
  combShannon$NormShan[combShannon$NormShan == 0] = minVal
  combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))
  
  
  # univariate nin parametric tests 
  
  #combShannon = combShannon[!combShannon$days_post_symptom_onset > 100,]
  combShannon = combShannon[!((combShannon$AP_lab_id == "EHC_C19_3662V") | (combShannon$AP_lab_id == "EHC_C19_3593E")),]
  
  combShannon = combShannon[combShannon$vax_doses_received < 3,]
  
  varList = c("Ct_depth_adj", "Ct_value","vax_doses_received", "meandepth", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
  varList = c("vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
  ShanCorAll = runSpearman(df=combShannon, varList=varList, testVar="Shannon", test = "spearman")
  wilcVax_2_All = runWilcox(curPred ="Vaccinated", combShannon= combShannon)
  wilcVax = runWilcox(curPred ="vax_doses_received", combShannon= combShannon)
  combShannon = combShannon[combShannon$WHO_variant != "Beta",]
  wilcVarAll = runWilcox(curPred ="WHO_variant", combShannon= combShannon)
  
  
  # multivariate model
  print(shapiro.test(combShannon$NormShan))
  print(ks.test(combShannon$NormShan, "pnorm"))
  
  DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
  varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
  varSel = varSel[!varSel$days_post_symptom_onset > 100,]
  varSel = drop_na(varSel)
  multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
  
  print(shapiro.test(varSel$NormShan))
  print(ks.test(varSel$NormShan, "pnorm"))
  
  # step AIC
  model1 = lm(NormShan~Protocol+meandepth+ Ct_depth_adj + WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) , data = varSel) 
  step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
  residuals <- residuals(step_model)
  print(ks.test(residuals, "pnorm"))
  print(shapiro.test(residuals))
  fullMult = summary(step_model)
  fullMult = data.frame(fullMult$coefficients)
  multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
  colnames(fullMult) = multColumns
  # write.csv(fullMult, "Paper/Tables/shannon_Multivar_Freq_1-0.02_addMiss.csv", row.names = T)
  
  print(fullMult)
  
  varSel = combShannon[combShannon$WHO_variant == "Delta" | combShannon$WHO_variant == "Omicron",]
  custMod = lm(NormShan~Protocol+Ct_depth_adj + WHO_variant+Vaccinated + days_post_symptom_onset , data = varSel) 
  custMod = lm(NormShan~Protocol+Ct_depth_adj + WHO_variant+Vaccinated , data = varSel) 
  custMult = summary(custMod)
  custMult = data.frame(custMult$coefficients)
  multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
  colnames(custMult) = multColumns
  
  print(custMult)
  
  if (addMiss == T & writeRes == T) {
    write.csv(ShanCorAll, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Spearman_Freq_1-",minFreq, "_addMiss_2023-12-20.csv"), row.names = F)
    write.csv(wilcVax_2_All, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Wilcox_Vax_Freq_1-",minFreq, "_addMiss_2023-12-20.csv"), row.names = F)
    write.csv(wilcVax, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Wilcox_VaxGroups_Freq_1-", minFreq, "_addMissv_2023-12-20.csv"), row.names = F)
    write.csv(wilcVarAll, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Wilcox_WHO_Freq_1-", minFreq, "_addMiss_2023-12-20.csv"),  row.names = F)
    
    write.csv(fullMult, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_StepMultivar_Freq_1-", minFreq, "_addMiss_2023-12-20.csv"), row.names = T)
    write.csv(custMult, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Multivar_Freq_1-", minFreq, "_addMiss_2023-12-20.csv"), row.names = T)
    
  } else if (addMiss == F & writeRes == T) {
    write.csv(ShanCorAll, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Spearman_Freq_1-",minFreq, "_2023-12-20.csv"), row.names = F)
    write.csv(wilcVax_2_All, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Wilcox_Vax_Freq_1-",minFreq, "_2023-12-20.csv"), row.names = F)
    write.csv(wilcVax, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Wilcox_VaxGroups_Freq_1-", minFreq, "_2023-12-20.csv"), row.names = F)
    write.csv(wilcVarAll, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Wilcox_WHO_Freq_1-", minFreq, "_2023-12-20.csv"),  row.names = F)
    
    write.csv(fullMult, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_StepMultivar_Freq_1-", minFreq, "_2023-12-20.csv"), row.names = T)
    write.csv(custMult, paste0("Paper/Tables/VaxDose_0-2/Shannon_MeanOverlap_GenLenAdj_Multivar_Freq_1-", minFreq, "_2023-12-20.csv"), row.names = T)
    
  }
  return(combShannon)
}

dir.create("Paper/Tables/VaxDose_0-2/")

#shanInd = makeAllTables(minFreq=0.01, addMiss=F, writeRes = T)
# shanInd = makeAllTables(minFreq=0, addMiss=F, writeRes = T)


#shanInd = makeAllTables(minFreq=0, addMiss=T, writeRes = T)
shanInd = makeAllTables(minFreq=0.01, addMiss=T, writeRes = F)


shanInd = shanInd[!is.na(shanInd$AP_lab_id),]

#write.csv(shanInd, "Paper/Tables/VaxDose_0-2/ShannInd_MeanOverlap_GenLenAdj_Freq_1-0.01_2023-12-20.csv", row.names = F)

system("aws s3 sync Paper/Tables/ s3://abombin/ARVAR/iSNVs/Paper/Tables/")

varSel = shanInd[shanInd$WHO_variant == "Delta" | shanInd$WHO_variant == "Omicron",]
custMod = lm(NormShan~Protocol+meandepth + Ct_value + WHO_variant+Vaccinated + days_post_symptom_onset , data = varSel) 
custMod = lm(NormShan~Protocol+meandepth + Ct_value  + WHO_variant+Vaccinated , data = varSel) 
custMult = summary(custMod)
