library(fitdistrplus)
library(nortest)
library(glmmTMB)
library(MASS)
library(rstatix)
library(dplyr)
library(pscl)
library(car)

dir.create("IntraSnv_results/statistic_res")

df = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv')
metadat = read.csv("Final_vaxbt_dataset_AP_metadata.csv")

df = df[df$Origin == "Overlap",]

metadatFilt = unique(metadat[, c("AP_lab_id", "Ct_value", "WHO_variant", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity" )])
metadatFilt[metadatFilt == ""] <- NA
metadatFilt[metadatFilt == "."] <- NA
#metadatFilt = metadatFilt[!is.na(metadatFilt$WHO_variant),]

metadatFilt$Sample1 = gsub('_', "-", metadatFilt$AP_lab_id)
metadatFilt$Sample = gsub('_', "-", metadatFilt$AP_lab_id)

combDat = plyr::join(df, metadatFilt, by = "Sample", type = "left", match = 'all')

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

combDat$Pi.Ln.Pi. = ShPi(targDf=combDat, freqCol="ALLELE.FREQUENCY")

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

shannonDf = getShannon(df=combDat)

combShannon = plyr::join(shannonDf, metadatFilt, by = "Sample", type = "left", match = 'all')

combShannon = combShannon[!is.na(combShannon$AP_lab_id),]


columns_to_check <- c("Ct_value", "vax_doses_received", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")
combShannon[columns_to_check] <- lapply(combShannon[columns_to_check], as.numeric)
columns_to_check <- c("Shannon", "meandepth", "Ct_value", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")


# Loop through the columns and perform Shapiro-Wilk test
for (col in columns_to_check) {
  p_value <- shapiro.test(combShannon[[col]])$p.value
  
  cat(paste("Shapiro-Wilk test for column", col, "\n"))
  cat("p-value:", p_value, "\n")
  
  if (p_value > 0.05) {
    cat("Conclusion: Data is approximately normally distributed (p > 0.05)\n\n")
  } else {
    cat("Conclusion: Data is not normally distributed (p <= 0.05)\n\n")
  }
}


var(combShannon$Shannon)
mean(combShannon$Shannon)

combShannon$NormShan = scale(log(combShannon$Shannon, base = 10))
combShannon$Ct_depth_adj = combShannon$Ct_value / combShannon$meandepth


model = glmmTMB(Shannon ~ Ct_value, data = combShannon, family = Gamma(link = "log"))
model <- glm.nb(Shannon ~ Ct_depth_adj, data = combShannon)
model <- glm.nb(Shannon ~ Ct_value, data = combShannon) 
model <- lm(NormShan~ Ct_value, data = combShannon) 
model <- lm(Shannon~ Ct_value, data = combShannon) #

summary(model)
# Obtain the residuals from the model
residuals <- residuals(model, type = "pearson")
# Perform a KS test on the residuals
ks.test(residuals, "pnorm", mean = 0, sd = 1)
shapiro.test(residuals)


# multivariate
model = glmmTMB(Shannon ~ meandepth + Ct_value + days_since_last_vax + days_post_symptom_onset + WHO_variant, data = combShannon, family = Gamma(link = "log"))
model <- glm.nb(Shannon ~ meandepth + Ct_value + days_since_last_vax + days_post_symptom_onset + WHO_variant, data = combShannon)
summary(model)
# Obtain the residuals from the model
residuals <- residuals(model, type = "pearson")
# Perform a KS test on the residuals
ks.test(residuals, "pnorm", mean = 0, sd = 1)
shapiro.test(residuals)
# 

runSpearman = function(df, varList, testVar) {
  Variable = character()
  Estimate = numeric()
  Pval = numeric()
  curResp = df[[testVar]]
  for (i in varList ) {
    #print(i)
    curPredict = df[[i]]
    curCor = cor.test(curResp, curPredict, method = "spearman")
    cur_cor_coef = curCor$estimate
    cur_pval = curCor$p.value
    Variable = c(Variable, i)
    Estimate = c(Estimate, cur_cor_coef)
    Pval = c(Pval, cur_pval)
  }
  combDat = data.frame(Variable, Estimate, Pval)
}

varList = c("Ct_depth_adj", "Ct_value","vax_doses_received", "meandepth", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")

shanCor = runSpearman(df=combShannon, varList=varList, testVar="Shannon")



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
  return(wilcox_res)
}

wilcVax = runWilcox(curPred ="vax_doses_received", combShannon= combShannon)
wilcVax = wilcVax[wilcVax$n1 >= 4 & wilcVax$n2 >= 4, ]
colnames(wilcVax)[1] = "Variable"
wilcVar = runWilcox(curPred ="WHO_variant", combShannon= combShannon)
wilcVar = wilcVar[wilcVar$n1 >= 4 & wilcVar$n2 >= 4, ]
colnames(wilcVar)[1] = "Variable"
sel_columns = c("Variable", "group1", "group2" , "p", "Groups_difference")

wilcVax=  wilcVax[, sel_columns]
wilcVar= wilcVar[, sel_columns]


# polymodel
model1  <- lm(NormShan ~ days_since_last_vax+I(days_since_last_vax^2), data = combShannon) # not sign
model1  <- lm(NormShan ~ meandepth+I(meandepth^2), data = combShannon) # sign
model1  <- lm(NormShan ~ days_post_symptom_onset+I(days_post_symptom_onset^2), data = combShannon) # sign partially norm
model1  <- lm(NormShan ~ vax_doses_received+I(vax_doses_received^2), data = combShannon) # sign
summary(model1)

residuals <- residuals(model1, type = "pearson")
# Perform a KS test on the residuals
ks.test(residuals, "pnorm", mean = 0, sd = 1)
shapiro.test(residuals)


# multivariate model
multModel = lm(NormShan ~ meandepth + Ct_value + vax_doses_received + days_since_last_vax + disease_severity, data = combShannon) # not normally distributed residuals
#multModel = lm(NormShan ~ meandepth + Ct_depth_adj + vax_doses_received + days_since_last_vax + disease_severity, data = combShannon)
#multModel = glm.nb(Shannon ~ Ct_depth_adj + vax_doses_received + days_since_last_vax + disease_severity, data = combShannon)
#multModel = glmmTMB(Shannon ~  Ct_depth_adj + vax_doses_received + days_since_last_vax + disease_severity, data = combShannon, family = Gamma(link = "log"))
summary(multModel)
# Obtain the residuals from the model
residuals <- residuals(multModel, type = "pearson")
# Perform a KS test on the residuals
ks.test(residuals, "pnorm", mean = 0, sd = 1)
shapiro.test(residuals)


stepAIC(multModel, direction = "both" , trace = T, steps = 10000)

model3  <- lm(NormShan ~ vax_doses_received, data = combShannon)
summary(model3)


## run Additional models for Anne

combShannon$Vaccinated = NA
combShannon$Vaccinated[combShannon$vax_doses_received == 0]<-"No"
combShannon$Vaccinated[combShannon$vax_doses_received > 0]<-"Yes"


wilcVax_2 = runWilcox(curPred ="Vaccinated", combShannon= combShannon)
colnames(wilcVax_2)[1] = "Variable"
wilcVax_2 =  wilcVax_2[, sel_columns]

combShannonSub = combShannon[combShannon$WHO_variant == "Delta" | combShannon$WHO_variant =="Omicron",]

multModel3 = lm(NormShan ~ Ct_value + meandepth + Vaccinated + days_post_symptom_onset + WHO_variant , data = combShannonSub)
summary(multModel3)
# Obtain the residuals from the model
residuals <- residuals(multModel3, type = "pearson")
# Perform a KS test on the residuals
ks.test(residuals, "pnorm", mean = 0, sd = 1)
shapiro.test(residuals)



multModel4 = lm(NormShan ~  Vaccinated + WHO_variant , data = combShannonSub)
summary(multModel4)
# Obtain the residuals from the model
residuals <- residuals(multModel4, type = "pearson")
# Perform a KS test on the residuals
ks.test(residuals, "pnorm", mean = 0, sd = 1)
shapiro.test(residuals)


fullMult = summary(multModel3)
fullMult = data.frame(fullMult$coefficients)
colnames(fullMult) = multColumns

partMult = summary(multModel4)
partMult = data.frame(partMult$coefficients)
colnames(partMult) = multColumns