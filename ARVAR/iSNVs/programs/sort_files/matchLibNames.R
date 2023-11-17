# pattern <- "EHC.*C19.*1412J.*L2.*"
# is_match <- grepl(q2, "EHC_C19-1412J-L2")

setwd("/home/ubuntu/extraVol/ARVAR/iSNVs")

metaseq_samples1 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/fastqs', pattern = "_R1")
metaseq_samples2 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-29/input', pattern = "_R1")

ampseq_samples1 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/input', pattern = "_R1")
ampseq_samples2 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/input', pattern = "_R1")

formatNames = function(samplesList) {
  newNames = character()
  sampleNames = gsub("-", "_", samplesList)
  for ( curSample in sampleNames ) {
    sampNamelist = strsplit(curSample, "_")
    sampleName =  sampNamelist[[1]][1:4]
    sampleName = paste(sampleName, collapse = "_")
    sampleName = gsub('_R1.fastq.gz', "", sampleName)
    sampleName = gsub('\\_S.*', "", sampleName)
    if (grepl("WATER", sampleName, ignore.case = T) == F) {
      newNames = c(newNames, sampleName)
    }
  }
  return(newNames)
}

metaseq_names1 = formatNames(metaseq_samples1)
metaseq_names2 = formatNames(metaseq_samples2)

comb_metaseq_names = unique(c(metaseq_names1, metaseq_names2))
length(comb_metaseq_names)

ampseq_names1 = formatNames(ampseq_samples1)
ampseq_names2 = formatNames(ampseq_samples2)

#ampseq_names2[!(ampseq_names2%in%ampseq_names1)]

comb_ampseq_names = unique(c(ampseq_names1, ampseq_names2))
length(comb_ampseq_names)

length(comb_ampseq_names[comb_ampseq_names%in%comb_metaseq_names])

myOverlap = comb_ampseq_names[comb_ampseq_names%in%comb_metaseq_names]

seqNamesDf = read.csv("SeqSamplesNames_Depth10.csv")
colnames(seqNamesDf)

addPresentSamples = function(df, metaList, ampList) {
  allMeta = character()
  allAmp = character()
  for ( i in 1:nrow(df) ) {
    curSample = df$Sample_id[i]
    curMeta = metaList[grepl(curSample, metaList)]
    if (length(curMeta) == 0) {
      curMeta = NA
    }
    curAmp = ampList[grepl(curSample, ampList)]
    if (length(curAmp) == 0) {
      curAmp = NA
    }
    
    metaString = paste(curMeta, collapse = ";")
    ampString = paste(curAmp, collapse = ";")
    
    allMeta = c(allMeta, metaString)
    allAmp = c(allAmp, ampString)
    
  }
  df$MetaFound = allMeta  
  df$AmpFound = allAmp
  return(df)
}
 
addMissing = function(df) {
  missingMeta = character()
  missingAmp = character()
  for ( i in 1:nrow(df) ) {
    targMeta = strsplit(df$MetaseqNames[i], ";")[[1]]
    curMeta = strsplit(df$MetaFound[i], ";")[[1]]
    curMisMeta = targMeta[!targMeta%in%curMeta]
    if (length(curMisMeta) == 0) {
      curMisMeta = NA
    }
    
    targAmp = strsplit(df$AmpseqNames[i], ";")[[1]]
    curAmp = strsplit(df$AmpFound[i], ";")[[1]]
    curMisAmp = targAmp[!targAmp%in%curAmp]
    if (length(curMisAmp) == 0) {
      curMisAmp = NA
    }
    
    missMetaStr = paste(curMisMeta, collapse = ";")
    missAmpStr = paste(curMisAmp, collapse = ";")
    
    missingMeta = c(missingMeta,  missMetaStr)
    missingAmp = c(missingAmp, missAmpStr)
  }
  df$Missing_Metaseq = missingMeta
  df$Missing_Ampseq = missingAmp
  return(df)
}

seqNamesFound = addPresentSamples(df = seqNamesDf, metaList = comb_metaseq_names, ampList = comb_ampseq_names)
colnames(seqNamesFound)

seqWMissing = addMissing(seqNamesFound)
colnames(seqWMissing)

formatDf = seqWMissing[, c("Sample_id", "MetaseqNames", "MetaFound", "Missing_Metaseq", "AmpseqNames", "AmpFound", "Missing_Ampseq")]
formatDf[formatDf == "NA"] <- NA

write.csv(formatDf, "SeqSamples_Eval_Depth10.csv", row.names = F)