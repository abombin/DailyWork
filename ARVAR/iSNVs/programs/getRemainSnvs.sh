cat incomplete_metaseq_all_samples.txt | parallel -j 8 'cd {} && python ../../programs/lofreq.py \
--allele_frequency 0 -r ../../programs/data/MN908947.3.fna \
-b ../../programs/bin -t 1'; 