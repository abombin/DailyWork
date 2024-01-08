
df = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions_oversamp.csv')

df1 = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions_oversamp.csv')


predPath = read.csv('IntraSnv_results/metaseq_comb_derep.csv')
outPath = read.csv('IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions.csv')

predPath = predPath[predPath$coverage >= 97,]

spurSNVs = predPath[!predPath$Samp_Pos_Ref_Alt%in%outPath$Samp_Pos_Ref_Alt,]

min(spurSNVs$ALLELE.FREQUENCY)
max(spurSNVs$ALLELE.FREQUENCY)
mean(spurSNVs$ALLELE.FREQUENCY)
median(spurSNVs$ALLELE.FREQUENCY)

print("----------------------")
min(outPath$ALLELE.FREQUENCY)
max(outPath$ALLELE.FREQUENCY)
mean(outPath$ALLELE.FREQUENCY)
median(outPath$ALLELE.FREQUENCY)

wilcox.test(outPath$ALLELE.FREQUENCY, spurSNVs$ALLELE.FREQUENCY)
t.test(outPath$ALLELE.FREQUENCY, spurSNVs$ALLELE.FREQUENCY)

# for ampseq filtered out SNVs are little bit lower in frequency
# for metaseq filtered out SNVs are little bit lower in frequency, even less than ampseq


# freq <= 0.01


predPath = read.csv('IntraSnv_results/ampseq_comb_derep.csv')
outPath = read.csv('IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv')

predPath = predPath[predPath$coverage >= 97,]
predPath = predPath[predPath$ALLELE.FREQUENCY >= 0.01,]
outPath = outPath[outPath$ALLELE.FREQUENCY >= 0.01,]

spurSNVs = predPath[!predPath$Samp_Pos_Ref_Alt%in%outPath$Samp_Pos_Ref_Alt,]

min(spurSNVs$ALLELE.FREQUENCY)
max(spurSNVs$ALLELE.FREQUENCY)
mean(spurSNVs$ALLELE.FREQUENCY)
median(spurSNVs$ALLELE.FREQUENCY)

print("----------------------")
min(outPath$ALLELE.FREQUENCY)
max(outPath$ALLELE.FREQUENCY)
mean(outPath$ALLELE.FREQUENCY)
median(outPath$ALLELE.FREQUENCY)

wilcox.test(outPath$ALLELE.FREQUENCY, spurSNVs$ALLELE.FREQUENCY)
t.test(outPath$ALLELE.FREQUENCY, spurSNVs$ALLELE.FREQUENCY)

# for ampseq filtered out SNVs are the same in  frequency
# for metaseq filtered out SNVs are little bit lower in frequency