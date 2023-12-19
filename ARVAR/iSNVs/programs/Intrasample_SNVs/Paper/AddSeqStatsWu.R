
dir.create("Paper/Additional_analysis/", recursive = T)

metaseq = read.csv("snvs_comb_res/metaseq_comb_derep_decont_v2.csv")
ampseq = read.csv("snvs_comb_res/ampseq_comb_derep_v2.csv")

metaseq$Samp_Pos_Ref_Alt = paste(metaseq$Sample1, metaseq$POSITION, metaseq$REF.NT, metaseq$VAR.NT, sep = "__" )
ampseq$Samp_Pos_Ref_Alt = paste(ampseq$Sample1, ampseq$POSITION, ampseq$REF.NT, ampseq$VAR.NT, sep = "__" )

comb_amp_cov = read.csv("Paper/ampseqConsSeqStats.csv")
#colnames(comb_amp_cov)[1] = "OrigName"
comb_amp_cov = unique(comb_amp_cov)

comb_meta_cov = read.csv("Paper/metaseqConsSeqStats.csv")
#colnames(comb_meta_cov)[1] = "OrigName"
comb_meta_cov = unique(comb_meta_cov)


metaseqCov = plyr::join(metaseq, comb_meta_cov, by="OrigName", type="left", match = "all")
any(is.na(metaseqCov$coverage))

ampseqCov = plyr::join(ampseq, comb_amp_cov, by="OrigName", type="left", match = "all")
any(is.na(ampseqCov$coverage))

ampseqCov = ampseqCov[ampseqCov$coverage >= 97,]
metaseqCov = metaseqCov[metaseqCov$coverage >= 97,]

ampseqCov$Sample = ampseqCov$Sample1
metaseqCov$Sample = metaseqCov$Sample1

write.csv(ampseqCov, "Paper/ampseqWu_derep_v3.csv")
write.csv(metaseqCov, "Paper/metaseqWu_derep_decont_v3.csv")

# get overlapping 

overlapSamples = unique(ampseqCov$Sample1[ampseqCov$Sample1%in%metaseqCov$Sample1])

ampseqOverlap = ampseqCov[ampseqCov$Sample1%in%overlapSamples,]
metaseqOverlap = metaseqCov[metaseqCov$Sample1%in%overlapSamples,]

getConsensus <- function(metaSeq, ampSeq, protocol, maxFreq, minFreq, freqCol, Snv_col) {
  metaseq <- metaSeq
  ampseq <- ampSeq
  metaseq =  metaseq[metaseq$coverage >= 97,]
  ampseq =   ampseq[ampseq$coverage >= 97,]
  metaseqFilt <- metaseq[metaseq[[freqCol]] >= minFreq & metaseq[[freqCol]] <= maxFreq, ]
  ampseqFilt <- ampseq[ampseq[[freqCol]] >= minFreq & ampseq[[freqCol]] <= maxFreq, ]
  ampseqFilt = ampseqFilt[ampseqFilt$Sample%in%metaseqFilt$Sample,]
  metaseqFilt = metaseqFilt[metaseqFilt$Sample%in%ampseqFilt$Sample,]
  if (protocol == "ampseq") {
    # dataframe with stats
    targDf = ampseqFilt
    # data to check agains
    refDf = metaseqFilt
    targSnv <- unique(refDf[, Snv_col])
  } else if (protocol == "metaseq") {
    # dataframe with stats
    targDf = metaseqFilt
    # data to check agains
    refDf = ampseqFilt
    targSnv <- unique(refDf[, Snv_col])
  }
  ConsTest <- numeric()
  for (i in 1:nrow(targDf)) {
    curSnv =  targDf[i, Snv_col]
    if (curSnv%in%targSnv){
      curCons = 1
    } else {
      curCons = 0
    }
    ConsTest = c(ConsTest, curCons)
  }
  targDf$ConsTest <- ConsTest
  return(targDf)
}

ampseqCons = getConsensus(metaSeq=metaseqOverlap, ampSeq=ampseqOverlap, protocol="ampseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY", Snv_col = "Samp_Pos_Ref_Alt")
metaseqCons = getConsensus(metaSeq=metaseqOverlap, ampSeq=ampseqOverlap, protocol="metaseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY", Snv_col = "Samp_Pos_Ref_Alt")

table(ampseqCons$ConsTest)
table(metaseqCons$ConsTest)

write.csv(ampseqCons ,"Paper/Additional_analysis/ampseqWu_ConsTest_freq_1_0.csv", row.names = F)
write.csv(metaseqCons ,"Paper/Additional_analysis/metaseqWu_ConsTest_freq_1_0.csv", row.names = F)