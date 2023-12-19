

# cd IntraSnv_ampseq_overlap; time ls -d */ | parallel -j 16 'cd {} && python ../../programs/src/getVarFilesCLT.py -i output.sam --mapqual 0 --basequal 0 --allele_frequency 0 \
# -r reference.fa -b ../../programs/bin -t 1'; cd ../; 


# cd IntraSnv_metaseq_overlap; time ls -d */ | parallel -j 24 'cd {} && python ../../programs/src/getVarFilesCLT.py -i output.sam --mapqual 0 --basequal 0 --allele_frequency 0 \
# -r reference.fa -b ../../programs/bin -t 1'; cd ../; 

cd IntraSnv_metaseq_all; time ls -d */ | parallel -j 29 'cd {} && python ../../programs/src/getVarFilesCLT.py -i output.sam --mapqual 0 --basequal 0 --allele_frequency 0 \
-r reference.fa -b ../../programs/bin -t 1'; cd ../; 