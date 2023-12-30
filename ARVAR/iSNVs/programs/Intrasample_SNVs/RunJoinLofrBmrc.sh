# cd IntraSnv_ampseq_overlap; time ls -d */ | parallel -j 10 'cd {} && python ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/JoinLofrBamread.py -i sample_lofreq-output.tsv \
#   -r sample_pos-filter.tsv'; cd ../; 
  
  
# cd IntraSnv_metaseq_overlap; time ls -d */ | parallel -j 60 'cd {} && python ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/JoinLofrBamread.py -i sample_lofreq-output.tsv \
#   -r sample_pos-filter.tsv'; cd ../; 


cd IntraSnv_ampseq_all; time ls -d */ | parallel -j 15 'cd {} && python ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/JoinLofrBamread.py -i sample_lofreq-output.tsv \
   -r sample_pos-filter.tsv'; cd ../; 
  
cd IntraSnv_metaseq_all; time ls -d */ | parallel -j 15 'cd {} && python ~/extraVol/ARVAR/iSNVs/programs/Intrasample_SNVs/JoinLofrBamread.py -i sample_lofreq-output.tsv \
  -r sample_pos-filter.tsv'; cd ../; 