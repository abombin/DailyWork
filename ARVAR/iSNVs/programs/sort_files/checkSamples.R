setwd("/home/ubuntu/extraVol/ARVAR/iSNVs")

metaseq_samples1 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/fastqs', pattern = "_R1")
metaseq_samples2 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-29/input', pattern = "_R1")

ampseq_samples1 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/input', pattern = "_R1")
ampseq_samples2 = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/input', pattern = "_R1")

ampseq_samples2[!(ampseq_samples2%in%ampseq_samples1)]


formatNames = function(samplesList) {
  newNames = character()
  sampleNames = toupper(gsub("_", "-", samplesList))
  for ( curSample in sampleNames ) {
    sampNamelist = strsplit(curSample, "-")
    sampleName =  sampNamelist[[1]][1:3]
    sampleName = paste(sampleName, collapse = "-")
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

ampseq_names2[!(ampseq_names2%in%ampseq_names1)]

comb_ampseq_names = unique(c(ampseq_names1, ampseq_names2))
length(comb_ampseq_names)

length(comb_ampseq_names[comb_ampseq_names%in%comb_metaseq_names])

myOverlap = comb_ampseq_names[comb_ampseq_names%in%comb_metaseq_names]

# compare with Ludy's list
ludy_df = read.csv("complete_samples_ludy.csv")
ludySamples = unique(ludy_df$AP_lab_id)
length(ludySamples)
ludyNames = unique(formatNames(ludySamples))
length(ludyNames)

length(ludyNames[ludyNames%in%comb_metaseq_names])
length(ludyNames[ludyNames%in%comb_ampseq_names])

allMySamples = unique(c(comb_metaseq_names, comb_ampseq_names))
length(allMySamples)
length(ludyNames[ludyNames%in%allMySamples])

#(length(ludyNames) - (length(ludyNames[ludyNames%in%comb_metaseq_names]) + length(ludyNames[ludyNames%in%comb_ampseq_names])))
missSampNumb = length(ludyNames) - length(ludyNames[ludyNames%in%allMySamples])
missSamples = unique(ludyNames[!(ludyNames%in%allMySamples)])
length(missSamples)
length(missSamples) == missSampNumb

missSampNumb / length(ludyNames) * 100

ludyNames[!(ludyNames%in%allMySamples)]

# 
presMetaseq = ludyNames[(ludyNames%in%comb_metaseq_names)]
presAmpseq = ludyNames[(ludyNames%in%comb_ampseq_names)]
presOverlap = ludyNames[(ludyNames%in%myOverlap)]

# save results
write.table(presMetaseq, "Andrei_metaseqSamples_2023-07-05.csv", row.names = F, col.names = F)
write.table(presAmpseq, "Andrei_ampseqSamples_2023-07-05.csv", row.names = F, col.names = F)
write.table(presOverlap, "Andrei_overlapSamples_2023-07-05.csv", row.names = F, col.names = F)
write.table(missSamples, "Andrei_missingSamples_2023-07-05.csv", row.names = F, col.names = F)

system("aws s3 cp Andrei_metaseqSamples_2023-07-05.csv s3://abombin/Ludy_vaxbt/")
system("aws s3 cp Andrei_ampseqSamples_2023-07-05.csv s3://abombin/Ludy_vaxbt/")
system("aws s3 cp Andrei_overlapSamples_2023-07-05.csv s3://abombin/Ludy_vaxbt/")
system("aws s3 cp Andrei_missingSamples_2023-07-05.csv s3://abombin/Ludy_vaxbt/")


# compare with Ludys sequenced samples
ludy_df = read.csv("ludy_sequenced_samples.csv")
ludySamples = unique(ludy_df$Sample_ID)
length(ludySamples)
ludyNames = unique(formatNames(ludySamples))
length(ludyNames)

length(ludyNames[ludyNames%in%comb_metaseq_names])
length(ludyNames[ludyNames%in%comb_ampseq_names])

allMySamples = unique(c(comb_metaseq_names, comb_ampseq_names))
length(allMySamples)
length(ludyNames[ludyNames%in%allMySamples])

#(length(ludyNames) - (length(ludyNames[ludyNames%in%comb_metaseq_names]) + length(ludyNames[ludyNames%in%comb_ampseq_names])))
missSampNumb = length(ludyNames) - length(ludyNames[ludyNames%in%allMySamples])
missSamples = unique(ludyNames[!(ludyNames%in%allMySamples)])
length(missSamples)
length(missSamples) == missSampNumb

missSampNumb / length(ludyNames) * 100

ludyNames[!(ludyNames%in%allMySamples)]

# 
presMetaseq = ludyNames[(ludyNames%in%comb_metaseq_names)]
presAmpseq = ludyNames[(ludyNames%in%comb_ampseq_names)]
presOverlap = ludyNames[(ludyNames%in%myOverlap)]

# filter ludy by depth
ludy_df = read.csv("ludy_sequenced_samples.csv")
ludy_df$Depth = as.numeric(as.character(ludy_df$Depth))
ludy_df = ludy_df[!is.na(ludy_df$Depth),]
ludy_df = ludy_df[ludy_df$Depth > 10 ,]
ludySamples = unique(ludy_df$Sample_ID)
length(ludySamples)
ludyNames = unique(formatNames(ludySamples))
length(ludyNames)

allMySamples = unique(c(comb_metaseq_names, comb_ampseq_names))
length(allMySamples)
length(ludyNames[ludyNames%in%allMySamples])

missSampNumb = length(ludyNames) - length(ludyNames[ludyNames%in%allMySamples])
missSamples = unique(ludyNames[!(ludyNames%in%allMySamples)])
length(missSamples)
length(missSamples) == missSampNumb