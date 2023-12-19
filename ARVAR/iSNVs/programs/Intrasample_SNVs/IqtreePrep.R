library(phylotools)

alignment = read.fasta("/home/ubuntu/extraVol/ARVAR/iSNVs/All_Samples_Mapping/All/iqtree/All_samples_align.fasta")

alignment = read.fasta("All_Samples_Mapping/All/All_samples_align.fasta")
curAlignment = alignment[!alignment$seq.name == "MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome",]

meta_vector <- rep("metaseq", times = 399)
amp_vector = rep("ampseq", times = 712)
combVect = c(amp_vector, meta_vector)

curAlignment$SeqType = combVect


curAlignment$seq.name = gsub(" MN908947.3", "", curAlignment$seq.name)

curAlignment$seq.name = paste0(curAlignment$seq.name, "_", curAlignment$SeqType)
curAlignment = curAlignment[, 1:2]

wuHo = alignment[1,]
wuHo$seq.name[1] = 'MN908947.3'

finalAlign = rbind(wuHo, curAlignment)

phylotools::dat2fasta(finalAlign, "/home/ubuntu/extraVol/ARVAR/iSNVs/All_Samples_Mapping/All/iqtree/All_samples_edit.align")
