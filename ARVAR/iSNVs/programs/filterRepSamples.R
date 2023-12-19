editNames = function(df) {
  combNames = character()
  for (i in 1:nrow(df)) {
    curName = df$OrigName[i]
    newName = gsub("_", "-", curName)
    newName = sub(".*EHC", "EHC", newName)
    nameList = strsplit(newName, "-")[[1]]
    nameList= nameList[1:3]
    newName  = paste(nameList , collapse = "-")
    combNames = c(combNames, newName)
  }
  df$Sample1 = combNames
  return(df)
}


amp_new = read.csv("snvs_comb_res/ampseq_found.csv")

amp_old = read.csv("snvs_comb_res/ampseq_old.csv")

cov_amp_new = read.csv("ampseqCovDepth.csv")
cov_amp_old = read.csv("ampseqOldSampCovDepth.csv")
comb_amp_cov = rbind(cov_amp_new, cov_amp_old)
colnames(comb_amp_cov)[1] = "OrigName"

colnames(amp_new)

ampCompComb = rbind(amp_new, amp_old)
ampCompComb = editNames(df=ampCompComb)

write.csv(ampCompComb, "snvs_comb_res/ampseq_all_v2.csv", row.names = F)

#amp_old_filt = amp_old[amp_old$ALLELE.FREQUENCY >= 0.02,]
#amp_new_filt = amp_new[amp_new$ALLELE.FREQUENCY >= 0.02,]

#ampComb = rbind(amp_new_filt,amp_old_filt)

libsDf = unique(ampCompComb[, c("OrigName", "ExactSample", "Sample1", "Batch" )])

sumSampFreq = data.frame(table(libsDf$Sample1))

repSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq > 1])
singSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq < 2])
length(repSamples)
length(singSamples)

#
selectMaxCovDepth = function(df, samplesList, sumStatDf) {
  combSamples = c()
  for (curSample in samplesList) {
    curDf = df[df$Sample1 == curSample,]
    curDfCov = plyr::join(curDf, sumStatDf, by = "OrigName", type = "left", match = "all")
    maxCov = max(curDfCov$Coverage)
    covDf = curDfCov[curDfCov$Coverage == maxCov,]
    if (nrow(covDf) > 1) {
      maxDepth = max(covDf$Mean_depth)
      depthDf = covDf[covDf$Mean_depth == maxDepth,]
      selSampleName = depthDf$OrigName[1]
    } else {
      selSampleName = covDf$OrigName[1]
    }
    combSamples = c(combSamples, selSampleName)
  }
  return(combSamples)
}

maxDepthSamples = selectMaxCovDepth(df=libsDf, samplesList=repSamples, sumStatDf = comb_amp_cov)

derepLib =  libsDf[libsDf$OrigName%in%maxDepthSamples,]
singLib = libsDf[libsDf$Sample1%in%singSamples,]
combSingLib = rbind(derepLib, singLib)

ampFilterLibs = ampCompComb[ampCompComb$OrigName%in%combSingLib$OrigName,]

q1 = unique(ampFilterLibs[, c("OrigName", "ExactSample", "Sample1", "Batch" )])
q2 = data.frame(table(q1$Sample1))

write.csv(ampFilterLibs, "snvs_comb_res/ampseq_comb_derep_v2.csv", row.names = F)

q3 = ampFilterLibs[ampFilterLibs$ALLELE.FREQUENCY >= 0.02 &  ampFilterLibs$ALLELE.FREQUENCY <= 0.98,]

## metaseq
contamList = read.table("metaseqContam.list")
contamList = contamList$V1
meta_new = read.csv("snvs_comb_res/metaseq_found.csv")

meta_old = read.csv("snvs_comb_res/metaseq_old.csv")

cov_meta_new = read.csv("metaseqCovDepth.csv")
cov_meta_old = read.csv("metaseqOldSampCovDepth.csv")
comb_meta_cov = rbind(cov_meta_new, cov_meta_old)
colnames(comb_meta_cov)[1] = "OrigName"
colnames(meta_new)
metaCompComb = rbind(meta_new, meta_old)
metaCompComb = editNames(df=metaCompComb)

write.csv(metaCompComb, "snvs_comb_res/metaseq_all_v2.csv", row.names = F)

libsDf = unique(metaCompComb[, c("OrigName", "ExactSample", "Sample1", "Batch" )])
# remove contaminated samples
remCont = function(contamList, libsDf) {
  exclSamples = character()
  contamList = gsub("_", ".", contamList)
  samplesList = libsDf$OrigName
  for (i in contamList) {
    curSamples = samplesList[grepl(i, samplesList)]
    exclSamples = c(exclSamples, curSamples)
  }
  return(exclSamples)
}

excludeSamples = remCont(contamList, libsDf)

libsDf = libsDf[!libsDf$OrigName%in%excludeSamples,]

sumSampFreq = data.frame(table(libsDf$Sample1))
repSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq > 1])
singSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq < 2])
length(repSamples)
length(singSamples)

maxDepthSamples = selectMaxCovDepth(df=libsDf, samplesList=repSamples, sumStatDf = comb_meta_cov)

derepLib =  libsDf[libsDf$OrigName%in%maxDepthSamples,]
singLib = libsDf[libsDf$Sample1%in%singSamples,]
combSingLib = rbind(derepLib, singLib)

metaFilterLibs = metaCompComb[metaCompComb$OrigName%in%combSingLib$OrigName,]
q1 = unique(metaFilterLibs[, c("OrigName", "ExactSample", "Sample1", "Batch" )])
q2 = data.frame(table(q1$Sample1))
write.csv(metaFilterLibs, "snvs_comb_res/metaseq_comb_derep_decont_v2.csv", row.names = F)

q3 = metaFilterLibs[metaFilterLibs$ALLELE.FREQUENCY >= 0.02 &  metaFilterLibs$ALLELE.FREQUENCY <= 0.98,]