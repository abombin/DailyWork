refDf = read.csv("test_consensus/iSNVs_final_Dara.csv")

df = refDf
colnames(refDf)
colName = "Allele.frequency.at.Day.0"
sampleName = "2884X"

getTrueSet = function(df, colName, sampleName) {
  refDfFilt = df[!is.na(df[colName]),]
  #refDfFilt = refDfFilt[refDfFilt[, colName] < 0.90,]
  refDfFilt$Samples = sampleName
  refDfFilt$Samp_Pos_Ref_Alt = paste(refDfFilt$Samples, refDfFilt$NT.POSITION, refDfFilt$REF.NT, refDfFilt$VAR.NT , sep = "__")
  snvs = unique(refDfFilt$Samp_Pos_Ref_Alt)
  return(snvs)
}

x = getTrueSet(df = refDf, colName =  "Allele.frequency.at.Day.0" , sampleName = "2884X")
y = getTrueSet(df = refDf, colName =  "Allele.Frequency.at.Day.31" , sampleName = "2885Y")
z = getTrueSet(df = refDf, colName =  "Allele.Frequency.at.Day.44" , sampleName = "2886Z")

allConsSnv = unique(c(x,y,z))
length(allConsSnv)

# get vivacity data
metaResDf = read.csv("test_consensus/vivacity_metaConsCheck.csv")

metaResDf$ConsTest = NA
metaResDf$ConsTest[!(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 0
metaResDf$ConsTest[(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 1

sumConsTest = data.frame(table(metaResDf$ConsTest))

colnames(metaResDf)

# models

# multivariate
multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + DEPTH + QUAL + Ref_Al_RelPos + Var_Al_Relpos, data = metaResDf, family = binomial)
summary(multivarModel)

multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + Ref_Al_RelPos + Var_Al_Relpos, data = metaResDf, family = binomial)
summary(multivarModel)


# multivariate
multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + DEPTH + QUAL + Position_test, data = metaResDf, family = binomial)
summary(multivarModel)

multivarRes = data.frame(summary(multivarModel)$coefficients)
colnames(multivarRes)[4] = "Pval"
multivarRes$Variable = rownames(multivarRes)

# run all univar models
vars = c('ALLELE.FREQUENCY', 'STRAND.BIAS',  'DEPTH', 'Ref_Al_RelPos', "QUAL", 'Var_Al_Relpos', "Position_test", "Freq_adj", "Type", "Pi.Ln.Pi.")

runUnivarMod = function(vars, df) {
  resultsDf  = data.frame(matrix(nrow = 0, ncol = 0))
  for (curVar in vars) {
    testVar = df[["ConsTest"]]
    indVar = df[[curVar]]
    univarModel = glm(testVar ~ indVar, family = binomial)
    sumUniv = summary(univarModel)
    sumUnivDf = data.frame(sumUniv$coefficients)
    dfSub =  sumUnivDf[2,]
    dfSub$Var = curVar
    resultsDf = rbind(resultsDf, dfSub)
  }
  return(resultsDf)
}

univarRes = runUnivarMod(vars = vars, df = metaResDf)
colnames(univarRes)[4] = "Pval"

#write.csv(univarRes, file = "test_consensus/vivacity_ManualCheck_univarRes.csv", row.names = F)

mean(metaResDf$STRAND.BIAS[metaResDf$ConsTest == 1])
median(metaResDf$STRAND.BIAS[metaResDf$ConsTest == 1])

write.csv(multivarRes, file = "test_consensus/vivacity_DaraSet_multivar.csv", row.names = F)