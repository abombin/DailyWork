library(mice)

combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")
# Choose the column you want to transform
combShannon$NormShan = combShannon$Shannon
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))
DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
varSel = drop_na(varSel)
multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
model2 = lm(NormShan~Protocol++WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) 
fullMult = summary(model2)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns
fullMult$bonferroni_p = p.adjust(fullMult$Pval, method = "bonferroni")
fullMult$fdr_p = p.adjust(fullMult$Pval, method = "fdr")
fullMult$Predictor = rownames(fullMult)
predictorNames = c("Intercept","Metagenomic Sequencing","Omicron Variant ","Vaccinated","Days post symptom onset","Days post symptom onset^2")
fullMult$Predictor = predictorNames
write.csv(fullMult, "Paper/Tables/Original_Multivar.csv", row.names = F)

# Imputing
combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")
colSums(is.na(combShannon))
completeDat = combShannon[, c("Sample", "Shannon", "meandepth", "AP_lab_id", "WHO_variant", "disease_severity", "Protocol")]
IncomplDat = combShannon[, c("Ct_value", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "Ct_depth_adj", "Vaccinated")]
mice_data <- mice(IncomplDat, m = 10, seed=33)
imputed_data <- complete(mice_data)
combShannon = cbind(completeDat, imputed_data)
combShannon$Vaccinated[combShannon$vax_doses_received == 0]<-"No"
combShannon$Vaccinated[combShannon$vax_doses_received > 0]<-"Yes"
combShannon$NormShan = combShannon$Shannon
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))

DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
varSel = drop_na(varSel)
multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
model2 = lm(NormShan~Protocol+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) 
fullMult = summary(model2)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns
fullMult$bonferroni_p = p.adjust(fullMult$Pval, method = "bonferroni")
fullMult$fdr_p = p.adjust(fullMult$Pval, method = "fdr")
fullMult$Predictor = rownames(fullMult)
predictorNames = c("Intercept","Metagenomic Sequencing","Omicron Variant ","Vaccinated","Days post symptom onset","Days post symptom onset^2")
fullMult$Predictor = predictorNames
write.csv(fullMult, "Paper/Tables/Imputed_Multivar_Shannon.csv", row.names = F)

system("aws s3 cp --recursive Paper/Tables/ s3://abombin/ARVAR/iSNVs/Paper/Tables/")