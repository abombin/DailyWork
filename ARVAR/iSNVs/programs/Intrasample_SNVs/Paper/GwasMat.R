
ampseq = read.csv("Paper/Additional_analysis/ampseq_Wu_ConsTest_freq_1_0_Predictions.csv")
metaseq = read.csv("Paper/Additional_analysis/metaseq_Wu_ConsTest_freq_1_0_Predictions.csv")

metaseqRef = read.csv("Paper/metaseqWu_derep_decont_v3.csv")
metaseqRef = metaseqRef[metaseqRef$coverage >= 97,]
ampseqRef = read.csv("Paper/ampseqWu_derep_v3.csv")
ampseqRef = ampseqRef[ampseqRef$coverage >=97,]

length(unique(ampseq$OrigName)) == length(unique(ampseqRef$OrigName))
length(unique(metaseq$OrigName)) == length(unique(metaseqRef$OrigName))



ampseq$Protocol = "Ampseq"
metaseq$Protocol = "Metaseq"

combDat = rbind(ampseq, metaseq)

combDat = combDat[combDat$ALLELE.FREQUENCY >= 0.01,]

combDatSel = combDat[, c("Sample", "OrigName", "Protocol", "Batch", "POSITION", "REF.NT", "VAR.NT")]

combDatSel$Samp_Prot = paste(combDatSel$Sample, combDatSel$Protocol, sep = "__")
combDatSel$OrigName_Prot = paste(combDatSel$OrigName, combDatSel$Protocol, sep = "__")
length(unique(combDatSel$Sample))
length(unique(combDatSel$Samp_Prot))
length(unique(combDatSel$OrigName_Prot))

# find duplicates
sampProt = unique(combDatSel[, c("Sample", "Protocol")])
sumFreq = data.frame(table(sampProt$Sample))
dulicateNames = sumFreq$Var1[sumFreq$Freq > 1]
singletonNames = sumFreq$Var1[sumFreq$Freq < 2]
length(dulicateNames) + length(singletonNames) == length(unique(sampProt$Sample))

# check that all duplicates indeed have 1 ampseq and 1 metaseq samples
dupDf = sampProt[sampProt$Sample%in%dulicateNames,]
table(dupDf$Protocol)
dupDf = dupDf[dupDf$Protocol == "Ampseq",]
dupNameProt = paste(dupDf$Sample, dupDf$Protocol, sep = "__")

# get deduplicates SNVVs
singletonSNVs = combDatSel[combDatSel$Sample%in%singletonNames,]
dedupSNVs = combDatSel[combDatSel$Samp_Prot%in%dupNameProt,]

combSNVs = rbind(singletonSNVs, dedupSNVs)

# review SNVs that have more than 1 variant per position
combSNVs$Pos_Var = paste(combSNVs$POSITION, combSNVs$VAR.NT, sep = "__")
combSNVs$Sample_Pos_Var = paste(combSNVs$Sample, combSNVs$POSITION, combSNVs$VAR.NT, sep = "__")
snvsPerSamp = data.frame(table(combSNVs$Sample_Pos_Var))
min(snvsPerSamp$Freq)
max(snvsPerSamp$Freq)
freqVar = snvsPerSamp[snvsPerSamp$Freq > 1,]
freqSampVars = combSNVs[combSNVs$Sample_Pos_Var%in%freqVar$Var1,]

# need to include reference allele to get all deletions correctly
combSNVs$Pos_Ref_Var = paste(combSNVs$POSITION, combSNVs$REF.NT, combSNVs$VAR.NT, sep = "__")
sumFreq = data.frame(table(combSNVs$Pos_Ref_Var))
min(sumFreq$Freq)
max(sumFreq$Freq)
mean(sumFreq$Freq)
median(sumFreq$Freq)


makeSNVsMatrix = function(df, maf, sumFreq) {
  df$Sample = gsub("-", "_", df$Sample)
  nSamples = length(unique(df$Sample))
  samplesList = unique(df$Sample)
  mafTr = maf * nSamples
  selSnvs = sumFreq$Var1[sumFreq$Freq >= mafTr]
  allSnvsDat = data.frame()
  for (curSnv in selSnvs) {
    snvDf = df[df$Pos_Ref_Var == curSnv,]
    curPosition = unique(snvDf$POSITION)
    curRef = unique(snvDf$REF.NT)
    curVar = unique(snvDf$VAR.NT)
    combSamples = character()
    combVals = numeric()
    for (curSample in samplesList) {
      if (curSample%in%snvDf$Sample) {
        curVal = 1 
      } else {
        curVal = 0
      }
      combSamples = c(combSamples, curSample) 
      combVals = c(combVals, curVal)
    }
    names(combVals) = combSamples
    curSnvMat = data.frame(t(combVals))
    curSnvMeta = data.frame(POSITION=curPosition, REF.NT= curRef, VAR.NT=curVar)
    curSnvCombDat = cbind(curSnvMeta, curSnvMat)
    allSnvsDat = rbind(allSnvsDat, curSnvCombDat)
  }
  return(allSnvsDat)
}

snvsMat = makeSNVsMatrix(df=combSNVs, maf = 0.01, sumFreq=sumFreq)

write.csv(snvsMat, "Paper/Additional_analysis/GWAS/SnvsMatMaf0.01_MinFreq0.01_Wu.csv", row.names = F)
write.csv(combSNVs, "Paper/Additional_analysis/GWAS/AmpMetaCombSnvsPostPredDedup_MinFreq0.01_Wu.csv", row.names = F)

# check that records are correct
snvsMat[1, 1:3]
sum(snvsMat[1, 4:ncol(snvsMat)])
snvsMat[2, 1:3]
sum(snvsMat[2, 4:ncol(snvsMat)])

snvsMat[200, 1:3]
sum(snvsMat[200, 4:ncol(snvsMat)])

sumFreq$Freq[sumFreq$Var1 == "10029__C__T"]
sumFreq$Freq[sumFreq$Var1 == "10046__G__GT"]
sumFreq$Freq[sumFreq$Var1 == "22995__C__A"]