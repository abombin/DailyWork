# cd ampseq_vivacity_found; ls -d * | parallel -j 80 'cd {} && python ../../programs/src/filterLofreqCLT.py -i sample_lofreq-output.tsv \
#   -r sample_pos-filter.tsv -a ../../'; cd ../
#   
# cd ampseq_vivacity_old; ls -d * | parallel -j 80 'cd {} && python ../../programs/src/filterLofreqCLT.py -i sample_lofreq-output.tsv \
#   -r sample_pos-filter.tsv -a ../../'; cd ../

cd metaseq_vivacity_found; ls -d * | parallel -j 25 'cd {} && python ../../programs/src/filterLofreqCLT.py -i sample_lofreq-output.tsv \
  -r sample_pos-filter.tsv -a ../../'; cd ../
  
cd metaseq_vivacity_old; ls -d * | parallel -j 25 'cd {} && python ../../programs/src/filterLofreqCLT.py -i sample_lofreq-output.tsv \
  -r sample_pos-filter.tsv -a ../../'; cd ../