combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")
colSums(is.na(combShannon))

MeanShann <- aggregate(cbind(Shannon, Ct_depth_adj) ~ Sample, data = combShannon , FUN = mean)

colnames(combShannon)

ref = unique(combShannon[, c("Sample", "WHO_variant", "vax_doses_received", "days_post_symptom_onset", "Vaccinated", "Protocol")])

JoinShannon = plyr::join(MeanShann, ref, by = "Sample", type = "left", match = "all")
combShannon=JoinShannon

# Choose the column you want to transform
combShannon$NormShan = combShannon$Shannon
min(combShannon$NormShan[combShannon$NormShan != 0])
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))

shapiro.test(combShannon$NormShan)
ks.test(combShannon$NormShan, "pnorm")


DfSel = combShannon[, c("Shannon", "NormShan", "WHO_variant", "vax_doses_received", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol")]
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated , data = varSel)
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)


varSel = drop_na(varSel)

model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) 
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

fullMult = summary(model2)

stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel )
residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)