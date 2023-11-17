library(mice)

combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")

colSums(is.na(combShannon))

completeDat = combShannon[, c("Sample", "Shannon", "meandepth", "AP_lab_id", "WHO_variant", "disease_severity", "Protocol")]
IncomplDat = combShannon[, c("Ct_value", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "Ct_depth_adj", "Vaccinated")]

mice_data <- mice(IncomplDat, m = 10, seed=33)

imputed_data <- complete(mice_data)

combShannon = cbind(completeDat, imputed_data)

combShannon$Vaccinated[combShannon$vax_doses_received == 0]<-"No"
combShannon$Vaccinated[combShannon$vax_doses_received > 0]<-"Yes"

colSums(is.na(combShannon))


combShannon$NormShan = combShannon$Shannon
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))
shapiro.test(combShannon$NormShan)
ks.test(combShannon$NormShan, "pnorm")


varSel = combShannon[combShannon$WHO_variant == "Delta" | combShannon$WHO_variant == "Omicron",]
nrow(varSel)
varSel = rstatix::drop_na(varSel)
nrow(varSel)

shapiro.test(varSel$NormShan)
ks.test(varSel$NormShan, "pnorm")

model2 = lm(NormShan~ Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel) # best fit but still not normal
model2 = lm(NormShan~ Protocol + Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel) 
model2 = lm(NormShan~ Protocol + Ct_depth_adj+WHO_variant+vax_doses_received + I(vax_doses_received^2) + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel) 
model2 = lm(NormShan~ Ct_depth_adj+WHO_variant+vax_doses_received + I(vax_doses_received^2) + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel) 


model2 = lm(NormShan~ Protocol + meandepth + I(meandepth^2) + Ct_value + I(Ct_value^2) +WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel) 
model2 = lm(NormShan~ Protocol + Ct_depth_adj + I(Ct_depth_adj^2) +WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel)


summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel )
residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
