metaseq = read.csv("snvs_comb_res/metaseq_comb_derep_decont_v2.csv")
ampseq = read.csv("snvs_comb_res/ampseq_comb_derep_v2.csv")

# 
# editNames = function(df) {
#   combNames = character()
#   for (i in 1:nrow(df)) {
#     curName = df$OrigName[i]
#     newName = gsub("_", "-", curName)
#     newName = sub(".*EHC", "EHC", newName)
#     nameList = strsplit(newName, "-")[[1]]
#     nameList= nameList[1:3]
#     newName  = paste(nameList , collapse = "-")
#     combNames = c(combNames, newName)
#   }
#   df$Sample1 = combNames
#   return(df)
# }
# 
# metaseq = editNames(df = metaseq)
# ampseq = editNames(df = ampseq)

metaseq$Samp_Pos_Ref_Alt = paste(metaseq$Sample1, metaseq$POSITION, metaseq$REF.NT, metaseq$VAR.NT, sep = "__" )
ampseq$Samp_Pos_Ref_Alt = paste(ampseq$Sample1, ampseq$POSITION, ampseq$REF.NT, ampseq$VAR.NT, sep = "__" )

cov_amp_new = read.csv("ampseqCovDepth.csv")
cov_amp_old = read.csv("ampseqOldSampCovDepth.csv")
comb_amp_cov = rbind(cov_amp_new, cov_amp_old)
colnames(comb_amp_cov)[1] = "OrigName"
comb_amp_cov = unique(comb_amp_cov)

cov_meta_new = read.csv("metaseqCovDepth.csv")
cov_meta_old = read.csv("metaseqOldSampCovDepth.csv")
comb_meta_cov = rbind(cov_meta_new, cov_meta_old)
colnames(comb_meta_cov)[1] = "OrigName"

metaseqCov = plyr::join(metaseq, comb_meta_cov, by="OrigName", type="left", match = "all")
any(is.na(metaseqCov$Coverage))

ampseqCov = plyr::join(ampseq, comb_amp_cov, by="OrigName", type="left", match = "all")
any(is.na(ampseqCov$Coverage))

metaseqCovFilt = metaseqCov[!metaseqCov$Coverage < 97,]
ampseqCovFit = ampseqCov[!ampseqCov$Coverage < 97,]

length(unique(metaseqCovFilt$Sample1))
length(unique(ampseqCovFit$Sample1))
length(unique(metaseqCovFilt$Sample1[metaseqCovFilt$Sample1%in%ampseqCovFit$Sample1]))
overlapSamples = unique(metaseqCovFilt$Sample1[metaseqCovFilt$Sample1%in%ampseqCovFit$Sample1])

write.csv(metaseqCovFilt, "snvs_comb_res/metaseq_comb_derep_decont_covFilt_97_v2.csv", row.names = F)
write.csv(ampseqCovFit, "snvs_comb_res/ampseq_comb_derep_covFilt_97_v2.csv", row.names = F)

# select only sample that overlap

metaOverlap = metaseqCovFilt[metaseqCovFilt$Sample1%in%overlapSamples,]
ampOverlap = ampseqCovFit[ampseqCovFit$Sample1%in%overlapSamples,]
length(unique(metaOverlap$Sample1))
length(unique(ampOverlap$Sample1))

write.csv(metaOverlap, "snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_97_v2.csv", row.names = F)
write.csv(ampOverlap, "snvs_comb_res/ampseq_overlap_comb_derep_covFilt_97_v2.csv", row.names = F)

# get overlapping set for all coverage
overlapSamplesAll = unique(metaseqCov$Sample1[metaseqCov$Sample1%in%ampseqCov$Sample1])

metaOverlapAll = metaseqCov[metaseqCov$Sample1%in%overlapSamplesAll,]
ampOverlapAll = ampseqCov[ampseqCov$Sample1%in%overlapSamplesAll,]

write.csv(metaOverlapAll, "snvs_comb_res/metaseq_overlap_comb_derep_decont_all_v2.csv", row.names = F)
write.csv(ampOverlapAll, "snvs_comb_res/ampseq_overlap_comb_derep_all_v2.csv", row.names = F)