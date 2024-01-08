# ampseq 

ampCompComb = read.csv("IntraSnv_results/ampseq_comb.csv")
comb_amp_cov = read.csv("IntraSnv_results/ampseq_stats_2.csv")
colnames(comb_amp_cov)[1] = "OrigName"

libsDf = unique(ampCompComb[, c("OrigName", "ExactSample", "Sample")])

sumSampFreq = data.frame(table(libsDf$Sample))

repSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq > 1])
singSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq < 2])
length(repSamples)
length(singSamples)

selectMaxCovDepth = function(df, samplesList, sumStatDf) {
  combSamples = c()
  for (curSample in samplesList) {
    curDf = df[df$Sample == curSample,]
    curDfCov = plyr::join(curDf, sumStatDf, by = "OrigName", type = "left", match = "all")
    maxCov = max(curDfCov$coverage)
    covDf = curDfCov[curDfCov$coverage == maxCov,]
    if (nrow(covDf) > 1) {
      maxDepth = max(covDf$meandepth)
      depthDf = covDf[covDf$meandepth == maxDepth,]
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
singLib = libsDf[libsDf$Sample%in%singSamples,]
combSingLib = rbind(derepLib, singLib)

ampFilterLibs = ampCompComb[ampCompComb$OrigName%in%combSingLib$OrigName,]

q1 = unique(ampFilterLibs[, c("OrigName", "ExactSample", "Sample" )])
q2 = data.frame(table(q1$Sample))

CombCov = plyr::join(ampFilterLibs, comb_amp_cov, by="OrigName", type="left", match = "all")
any(is.na(CombCov$coverage))

nrow(ampFilterLibs) == nrow(CombCov)

write.csv(CombCov, "IntraSnv_results/ampseq_comb_derep.csv", row.names = F)


#metaseq 

ampCompComb = read.csv("IntraSnv_results/metaseq_comb.csv")
comb_amp_cov = read.csv("IntraSnv_results/metaseq_stats_2.csv")
colnames(comb_amp_cov)[1] = "OrigName"

libsDf = unique(ampCompComb[, c("OrigName", "ExactSample", "Sample")])

sumSampFreq = data.frame(table(libsDf$Sample))

repSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq > 1])
singSamples = as.character(sumSampFreq$Var1[sumSampFreq$Freq < 2])
length(repSamples)
length(singSamples)

maxDepthSamples = selectMaxCovDepth(df=libsDf, samplesList=repSamples, sumStatDf = comb_amp_cov)

derepLib =  libsDf[libsDf$OrigName%in%maxDepthSamples,]
singLib = libsDf[libsDf$Sample%in%singSamples,]
combSingLib = rbind(derepLib, singLib)

ampFilterLibs = ampCompComb[ampCompComb$OrigName%in%combSingLib$OrigName,]

q1 = unique(ampFilterLibs[, c("OrigName", "ExactSample", "Sample" )])
q2 = data.frame(table(q1$Sample))

CombCov = plyr::join(ampFilterLibs, comb_amp_cov, by="OrigName", type="left", match = "all")
any(is.na(CombCov$coverage))

nrow(ampFilterLibs) == nrow(CombCov)

write.csv(CombCov, "IntraSnv_results/metaseq_comb_derep.csv", row.names = F)


# write only overlapping samples
rm(list=ls())
ampseq = read.csv('IntraSnv_results/ampseq_comb_derep.csv')
metaseq = read.csv("IntraSnv_results/metaseq_comb_derep.csv")

ampseqOverlap = ampseq[ampseq$Sample %in% metaseq$Sample,]
length(unique(ampseqOverlap$Sample))

metaseqOverlap = metaseq[metaseq$Sample %in% ampseq$Sample,]
length(unique(metaseqOverlap$Sample))

write.csv(ampseqOverlap, "IntraSnv_results/ampseq_comb_derep_overlap.csv", row.names = F)
write.csv(metaseqOverlap, "IntraSnv_results/metaseq_comb_derep_overlap.csv", row.names = F)