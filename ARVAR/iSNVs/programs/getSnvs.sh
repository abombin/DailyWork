
# newly found ampseq
# cd ampseq_vivacity_found; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/bamReadCount.py \
# --mapqual 0 --basequal 0 --allele_frequency 0 \
# -r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../;
#
# cd ampseq_vivacity_found; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/lofreq.py \
#  --allele_frequency 0 -r ../../programs/data/MN908947.3.fna \
#  -b ../../programs/bin -t 1'; cd ../



# previously found ampseq

# cd ampseq_vivacity_old; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/bamReadCount.py \
# --mapqual 0 --basequal 0 --allele_frequency 0 \
# -r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../; 
# 
# cd ampseq_vivacity_old; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/lofreq.py \
# --allele_frequency 0 -r ../../programs/data/MN908947.3.fna \
# -b ../../programs/bin -t 1'; cd ../
# 
# cd metaseq_vivacity_old; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/bamReadCount.py \
# --mapqual 0 --basequal 0 --allele_frequency 0 \
# -r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../; 
# 
# cd metaseq_vivacity_old; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/lofreq.py \
# --allele_frequency 0 -r ../../programs/data/MN908947.3.fna \
# -b ../../programs/bin -t 1'; cd ../

cd metaseq_vivacity_found; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/bamReadCount.py \
--mapqual 0 --basequal 0 --allele_frequency 0 \
-r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../;

cd metaseq_vivacity_found; time ls -d */ | parallel -j 94 'cd {} && python ../../programs/lofreq.py \
--allele_frequency 0 -r ../../programs/data/MN908947.3.fna \
-b ../../programs/bin -t 1'; cd ../

