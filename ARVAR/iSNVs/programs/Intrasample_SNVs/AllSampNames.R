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

colnames(amp_new)

ampCompComb = rbind(amp_new, amp_old)
ampCompComb = editNames(df=ampCompComb)
ampCombSamples = unique(ampCompComb[, c("Sample1", "OrigName")])

write.csv(ampCompComb, "snvs_comb_res/ampseq_all_v2.csv", row.names = F)
write.csv(ampCombSamples, "snvs_comb_res/ampseq_all_v2_names.csv", row.names = F)

## metaseq
meta_new = read.csv("snvs_comb_res/metaseq_found.csv")
meta_old = read.csv("snvs_comb_res/metaseq_old.csv")

metaCompComb = rbind(meta_new, meta_old)
metaCompComb = editNames(df=metaCompComb)

metaCombSamples = unique(metaCompComb[, c("Sample1", "OrigName")])

write.csv(metaCompComb, "snvs_comb_res/metaseq_all_v2.csv", row.names = F)
write.csv(metaCombSamples, "snvs_comb_res/metaseq_all_v2_names.csv", row.names = F)