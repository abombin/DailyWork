library(rstatix)
library(MASS)
library(car)
library(lme4)

combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")

combShannon$ShannYJ = yjPower(combShannon$Shannon, lambda = 0.5)

shapiro.test(combShannon$ShannYJ)
ks.test(combShannon$ShannYJ, "pnorm")


combShanFilt = combShannon[combShannon$Shannon > 0 ,]

# Choose the column you want to transform
combShanFilt$NormShan = scale(log(combShanFilt$Shannon, base = 10))

combShanFilt$NormShan = log(combShanFilt$Shannon, base = 10)

shapiro.test(combShanFilt$NormShan)
ks.test(combShanFilt$NormShan, "pnorm")

min(combShanFilt$Shannon) * 0.1

# min max normalization 

shannon <- combShannon$Shannon

# Min-max normalization
min_value <- min(shannon)
max_value <- max(shannon)

normalized_shannon <- (shannon - min_value) / (max_value - min_value)

# Replace the original Shannon column with the normalized values
combShannon$MinMaxShannon <- normalized_shannon
combShannon$MinMaxShannon = as.numeric(combShannon$MinMaxShannon)

shapiro.test(combShannon$MinMaxShannon)
ks.test(combShannon$MinMaxShannon, "pnorm")


###
colSums(is.na(combShannon))

combShannon$NormShan = combShannon$Shannon

combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
min(combShannon$Shannon)

combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))
shapiro.test(combShannon$NormShan)
ks.test(combShannon$NormShan, "pnorm")

# try models
colnames(combShannon)

model1 = lm(NormShan~ meandepth + Ct_value + days_post_symptom_onset, data = combShannon)
model1 = lm(MinMaxShannon~ meandepth + Ct_value + days_post_symptom_onset, data = combShannon)
summary(model1)
residuals <- residuals(model1)
ks.test(residuals, "pnorm")
shapiro.test(residuals)


# select variables with least NAs

DfSel = combShannon[, c("Shannon", "NormShan", "MinMaxShannon", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated")]
DfSel = drop_na(DfSel)

model1 = lm(NormShan~meandepth+Ct_value+WHO_variant+Vaccinated+disease_severity, data = DfSel)
summary(model1)
stepModel = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(stepModel)

residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

model2 = lm(NormShan~+Ct_value+WHO_variant+Vaccinated, data = DfSel) # good results but not normally distributed
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

# try only with Delta and Omicron
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
model2 = lm(NormShan~+Ct_value+WHO_variant+Vaccinated, data = varSel) # good results but not normally distributed
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
# try another combination of Ct and depth evaluation
model1 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated+disease_severity, data = DfSel)
summary(model1)
stepModel = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(stepModel)

residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

# check only with delta and omicron

varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
model2 = lm(NormShan~+Ct_depth_adj+WHO_variant+Vaccinated, data = varSel) # good results but not normally distributed
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepAIC(model2, direction = "both" , trace = T, steps = 10000)

##
#model3 <- glmer(Shannon ~ Ct_depth_adj+WHO_variant+Vaccinated, data = varSel, family = Gamma(link = log))
#model3 <- glm(Shannon  ~ Ct_depth_adj+WHO_variant+Vaccinated, family = poisson(link = "log"), data = varSel) # bad
varSel$Shannon[varSel$Shannon == 0] = 0.005900276
varSel$MinMaxShannon[varSel$MinMaxShannon < 0.000001] = 0.005900276
model3 = glmmTMB(MinMaxShannon ~  Ct_depth_adj + WHO_variant+Vaccinated, data = varSel, family = Gamma(link = "log"))
#model3 = glm(Shannon ~ Ct_depth_adj+WHO_variant+Vaccinated, family = Gamma(link = "log"), data = varSel) bad
residuals <- residuals(model3)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
summary(model3)

# try to renormalize only with selected groups
DfSel = combShannon[, c("Shannon", "NormShan", "MinMaxShannon", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated")]
DfSel = drop_na(DfSel)
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]

varSel$NormShan = varSel$Shannon
varSel$NormShan[varSel$NormShan == 0] = 0.005900276
min(varSel$NormShan)
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))
shapiro.test(combShannon$NormShan)
ks.test(combShannon$NormShan, "pnorm")

model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated, data = varSel) # results are worse than original
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepAIC(model2, direction = "both" , trace = T, steps = 10000)


model2 = lm(NormShan~Ct_depth_adj+WHO_variant+vax_doses_received + I(vax_doses_received^2), data = varSel) # results are worse than original
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepAIC(model2, direction = "both" , trace = T, steps = 10000)


# try with complete dataset 
varSel_all = combShannon[combShannon$WHO_variant == "Delta" | combShannon$WHO_variant == "Omicron",]
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated, data = varSel_all) # good results but not normally distributed
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepAIC(model2, direction = "both" , trace = T, steps = 10000)

# try to add more variables
colSums(is.na(combShannon))

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
stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel )
residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)


varSel = drop_na(varSel)
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) # ok results and somewhat normally distributed
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+vax_doses_received + I(vax_doses_received^2) + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) # the sameok results and somewhat normally distributed
model2 = lm(NormShan~Ct_depth_adj+WHO_variant+vax_doses_received + I(vax_doses_received^2) + days_post_symptom_onset , data = varSel) # better fit but more controversial
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel )
residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

model2 = lm(NormShan~Protocol +Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel) # protocol makes fit worse
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel)
residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

## test disease severity
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
model2 = lm(disease_severity~ NormShan + Ct_depth_adj+WHO_variant+Vaccinated , data = varSel)
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

varSel = drop_na(varSel)
model2 = lm(disease_severity~ NormShan+Ct_depth_adj+WHO_variant+Vaccinated + days_post_symptom_onset + I(days_post_symptom_onset^2), data = varSel)
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel)
residuals <- residuals(stepModel)
ks.test(residuals, "pnorm")
shapiro.test(residuals)

# test complete model
combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")
colSums(is.na(combShannon))
# Choose the column you want to transform
combShannon$NormShan = combShannon$Shannon
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))
DfSel = combShannon[, c("Shannon", "NormShan", "meandepth", "Ct_value", "WHO_variant", "vax_doses_received", "disease_severity", "Ct_depth_adj", "Vaccinated", "days_post_symptom_onset", "Protocol", "days_since_last_vax")]
#DfSel = drop_na(DfSel)
varSel = DfSel[DfSel$WHO_variant == "Delta" | DfSel$WHO_variant == "Omicron",]
varSel = drop_na(varSel)
model2 = lm(NormShan~Ct_depth_adj+WHO_variant + vax_doses_received + I(vax_doses_received^2) + days_post_symptom_onset + I(days_post_symptom_onset^2) + days_since_last_vax + I(days_since_last_vax^2), data = varSel)
summary(model2)
residuals <- residuals(model2)
ks.test(residuals, "pnorm")
shapiro.test(residuals)
stepModel = stepAIC(model2, direction = "both" , trace = T, steps = 10000)
summary(stepModel)
