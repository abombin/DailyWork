cd bactopia_gtdbtk
for i in $(cat ../bacteria.list); do
cd ./"$i"
mkdir -p select_output

for sample in $(cat samples.list); do

mkdir -p select_output/"$sample"/summary
mkdir -p select_output/"$sample"/antimicrobial-resistance
mkdir -p select_output/"$sample"/contigs
mkdir -p select_output/"$sample"/annotation
mkdir -p select_output/"$sample"/mapRefGen

cp ./output/"$sample"/quality-control/summary/*.html select_output/"$sample"/summary/
cp ./output/"$sample"/quality-control/summary/*.json select_output/"$sample"/summary/

cp ./output/"$sample"/antimicrobial-resistance/*  select_output/"$sample"/antimicrobial-resistance

cp ./output/"$sample"/assembly/*.fna.gz select_output/"$sample"/contigs/

cp ./output/"$sample"/annotation/*.tsv select_output/"$sample"/annotation
cp ./output/"$sample"/annotation/*.txt select_output/"$sample"/annotation

cp -r ./output/"$sample"/variants/ select_output/"$sample"/mapRefGen
cp -r ./output/"$sample"/variants/ select_output/"$sample"/mapRefGen
done
cd ../
done

echo "Bactopia Output Coppied"
sleep 20s

for i in $(cat ../bacteria.list); do
cd ./"$i"
mkdir -p select_output/microreact
cp ./panOut/bactopia-tools/pangenome/pangenome/iqtree/core-genome.contree ./select_output/microreact/"$i"_consensus.nwk
awk -F, 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' samples.list ../../metaData/metaDat.csv > metaDat.csv
cd ../
done

cd ../

echo "Coppied data for microreact"
sleep 20s


# compress output files tor transfer
now="$(date +'%m-%d-%Y')"
for i in $(cat ./bacteria.list)
do 
	#bucket=$(cat bacteria.bucket)
	zip -r ./bactopia_gtdbtk/"$i"/"$i"_"$now".zip ./bactopia_gtdbtk/"$i"/select_output
done



