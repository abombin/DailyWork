for i in $(cat ./sample.name)
do 
#/home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 30 --length 36 --paired "$i"_1.fastq.gz "$i"_2.fastq.gz
/home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20  --paired "$i"_R1.fastq.gz "$i"_R2.fastq.gz
done
