combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")
colSums(is.na(combShannon))

# Choose the column you want to transform
combShannon$NormShan = combShannon$Shannon
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))

shapiro.test(combShannon$NormShan)
ks.test(combShannon$NormShan, "pnorm")

DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
#DfSel = drop_na(DfSel)
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated , data = varSel)
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

multColumns = c("Estimate", "Std.Error", 't.value', "Pval")
fullMult = summary(model2)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns
fullMult$bonferroni_p = p.adjust(fullMult$Pval, method = "bonferroni")
fullMult$fdr_p = p.adjust(fullMult$Pval, method = "fdr")
fullMult$Predictor = rownames(fullMult)

#write.csv(fullMult, "IntraSnv_results/statistic_res/AllSamplesCombShannon_Multivar_CtVarVac.csv", row.names = F)

DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
varSel = drop_na(varSel)
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) 
fullMult = summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns
fullMult$bonferroni_p = p.adjust(fullMult$Pval, method = "bonferroni")
fullMult$fdr_p = p.adjust(fullMult$Pval, method = "fdr")
fullMult$Predictor = rownames(fullMult)

write.csv(fullMult, "IntraSnv_results/statistic_res/AllSamplesCombShannon_Multivar_CtVarVacSympt.csv", row.names = F)


model2 = lm(NormShan~Ct_depth_adj+WHO_variant+vax_doses_received + I(vax_doses_received^2) + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) 
fullMult = summary(model2)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns
fullMult$bonferroni_p = p.adjust(fullMult$Pval, method = "bonferroni")
fullMult$fdr_p = p.adjust(fullMult$Pval, method = "fdr")
fullMult$Predictor = rownames(fullMult)

write.csv(fullMult, "IntraSnv_results/statistic_res/AllSamplesCombShannon_Multivar_CtVarVacContSympt.csv", row.names = F)

system("aws s3 cp --recursive IntraSnv_results/statistic_res/ s3://abombin/ARVAR/iSNVs/IntraSnv_results/statistic_res/")