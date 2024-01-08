ampseq = read.csv('IntraSnv_results/ampseq_comb_derep.csv')
metaseq = read.csv('IntraSnv_results/metaseq_comb_derep.csv')

#ampseq = read.csv("/home/ubuntu/extraVol/ARVAR/iSNVs/IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_GlobPos.csv")
#metaseq = read.csv("/home/ubuntu/extraVol/ARVAR/iSNVs/IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_GlobPos.csv")

ampseq = ampseq[ampseq$coverage >= 97,]
metaseq = metaseq[metaseq$coverage >= 97, ]

length(unique(ampseq$OrigName))
length(unique(metaseq$OrigName))

q1 = unique(ampseq$OrigName)
q2 = unique(ampseq_ConsTest_freq_1_0$OrigName)

length(q1[!q1%in%q2])


packageVersion("stats")
packageVersion("MASS")
packageVersion("car")
packageVersion("mice")
packageVersion("ggplot2")

citation("MASS")
citation("car")
citation("mice")
citation("ggplot2")

# sum stats
length(unique(ampseq$Samp_Pos_Ref_Alt))
length(unique(metaseq$Samp_Pos_Ref_Alt))

ampseqFreq = data.frame(table(ampseq$OrigName))
mean(ampseqFreq$Freq)
median(ampseqFreq$Freq)
metaseqFreq = data.frame(table(metaseq$OrigName))
mean(metaseqFreq$Freq)
median(metaseqFreq$Freq)

nrow(ampseq[ampseq$ALLELE.FREQUENCY < 0.01 ,]) / length(unique(ampseq$Samp_Pos_Ref_Alt)) * 100
nrow(metaseq[metaseq$ALLELE.FREQUENCY < 0.01 ,]) / length(unique(metaseq$Samp_Pos_Ref_Alt)) * 100

# overlapping stats 
ampseq = read.csv("IntraSnv_results/ampseq_ConsTest_freq_1_0.csv")
metaseq = read.csv("IntraSnv_results/metaseq_ConsTest_freq_1_0.csv")

length(unique(ampseq$Sample_AlignPos_Ref_Var))
length(unique(metaseq$Sample_AlignPos_Ref_Var))

length(unique(ampseq$Sample_AlignPos_Ref_Var[ampseq$Sample_AlignPos_Ref_Var%in%metaseq$Sample_AlignPos_Ref_Var]))

# Random forest filtered SNVs summaries
metaseq = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions.csv')
ampseq = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv')


length(unique(ampseq$Samp_Pos_Ref_Alt))
length(unique(metaseq$Samp_Pos_Ref_Alt))

ampseqFreq = data.frame(table(ampseq$OrigName))
mean(ampseqFreq$Freq)
median(ampseqFreq$Freq)
metaseqFreq = data.frame(table(metaseq$OrigName))
mean(metaseqFreq$Freq)
median(metaseqFreq$Freq)
