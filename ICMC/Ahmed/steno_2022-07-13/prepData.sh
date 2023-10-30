mkdir -p bactopia_gtdbtk
mkdir -p test_gtdbtk/buckets
mkdir -p test_gtdbtk/taxa
mkdir -p input
mkdir -p original_files/R1
mkdir -p original_files/R2
mkdir -p gtdbtk/input

echo " Renaming files"

for i in ./original_files/*_R1_001.fastq.gz; do mv "$i" ./original_files/R1/; done
for i in ./original_files/*_R2_001.fastq.gz; do mv "$i" ./original_files/R2/; done
(cd ./original_files/R1 && ls *_R1_001.fastq.gz) > inputList.txt
#cat inputList.txt |cut -f1 -d"_" > newdir.txt # remove everything after first _
cp inputList.txt ./newdir.txt
sed -i 's/\_R1_001.fastq.gz//g' newdir.txt # remove specific string

 
#(cd ./input/ && ls *_R1_001.fastq.gz) > inputList.txt



for i in $(cat newdir.txt); do cp ./original_files/R1/"$i"_* ./input/"$i"_R1.fastq.gz; done
for i in $(cat newdir.txt); do cp ./original_files/R2/"$i"_* ./input/"$i"_R2.fastq.gz; done

#sed -i.bak 's/_R1.*//g' inputList.txt
#mv inputList.txt newdir.txt




for i in $(cat newdir.txt); do mkdir -p process_par/"$i"; mv ./input/"$i"_* ./process_par/"$i"/
  echo "$i" > process_par/"$i"/sample.name; done

echo "Prepared all files for pipeline"
sleep 10s