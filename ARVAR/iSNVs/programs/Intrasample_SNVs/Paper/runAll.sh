# complete pipeline to do analysis for the paper

# add EBS
sh start.sh

# go to directory
cd /home/ubuntu/extraVol/ARVAR/iSNVs

## Preprocessing 
# prepare all directories for metaseq and make alignments against the consensus genomes
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/AlignMetaseqAll.R
# prepare all directories for ampseq and make alignments against the consensus genomes
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/AlignOverlapAmpseqAll.R
# get SNV tables and annotations
sh ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/getSnvs.sh
# merge lofreq  and bamreadcount results
sh ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/RunJoinLofrBmrc.sh
# combine SNV tables for individual samples together in one table
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/combSnvsAll.R
# get sequencing stats for samples
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/sumDepthCov.R
# annotate SNVs with sequencing stats and remove duplicate samples per sequencing method
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/FilterRepSamples.R
# make alignments to map positions between overlapping samples
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/MapAllSamplesAlignment2.R
# Add alignment positions to overlapping samples
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/AddAlignPosCons2.R
## Filtering SNVs and Shannon index
# test default random forest prediction model for overlapping samples and make tables and AUC plot
python /home/ubuntu/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/RandomForest.py
# make predictions and filter the whole dataset with a default random forest method
python /home/ubuntu/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/RandForestRealPredict.py
# make alignment for all samples that passed random forest filtering to add SNVs positions between samples
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/MultiseqAlign.R
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/MultiseqIndicies.R
# add all samples alignment positions to the results of SNVs filtering
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/AddAllSamplesGlobPos.R
# Calculate Shannon index and make univariate and multivariate analysis, currently runs the datasets for oversampling method, so adjust input and output paths for default random forest
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/ShanTestsFreqFiltAverOverlapGenAdjImpute.R
## GWAS
# for the GWAS we use SNVs called against Wuhan reference
# get SNVs
sh ~/extraVol/ARVAR/iSNVs/programs/getSnvs.sh
# cpmbine lofreq and bamreadcount outputs and add SNVs annotatioins
sh ~/extraVol/ARVAR/iSNVs/programs/annotatePar.sh
# combine SNVs from individual samples to one table
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/combSnvDat.R
# get depth and coverage for consensus sequences, inefficient but curently needed to run the pipeline
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/sumDepthCov.R
# get all sequencing statistics for consensus sequences
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/ConsSeqStats.R
# remove duplicate samples per sequencing method
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/filterRepSamples.R
# remove samples with coverage less than 97%
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/filterCov.R
# add all sequencing stats to SNVs data
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/AddSeqStatsWu.R
# Test standard random forest model for filtering SNVs
python ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/RandForWu.py
# make real predictions with for all SNVs
python ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/RandForRealPRedWu.py
# copy ampseq and metaseq consensus sequences
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/CpConsMeta.R
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/CpConsAmp.R
# the following programs are set for oversampling random forest filtering strategy, change for default strategy in future
# make a data matrix with SNVs for GWAS
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/GwasMat.R
# make alignment of the consensus sequences to calculate the distances, and iqtree that is curently not used in the pipeline
sh ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/AlignIqtreeConsWu.sh
# calculate mds for population adjustment 
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/Pheno_Covar.R
# make GWAS with binary vaccination status
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/GwasVax.R
# GWAS with vaccination status as continius variable
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/GwasVaxDoses.R
# GWAS with disease severity 
Rscript --vanilla ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/Paper/GwasDis.R





