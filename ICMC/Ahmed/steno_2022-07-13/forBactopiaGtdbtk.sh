(cd ./test_gtdbtk/buckets && ls ) > ./test_gtdbtk/bacteria.list


for bucket in $(cat ./test_gtdbtk/bacteria.list)
do
mkdir -p ./bactopia_gtdbtk/"$bucket"/input
for sample in $(cat ./test_gtdbtk/buckets/"$bucket")
do
mv ./process_par/"$sample"/*.fastq.gz ./bactopia_gtdbtk/"$bucket"/input
echo "$sample" >> ./bactopia_gtdbtk/"$bucket"/samples.list
done 
done

echo "samples moved"
echo ""

for bucket in $(cat ./test_gtdbtk/bacteria.list)
do
cp ./test_gtdbtk/taxa/"$bucket" ./bactopia_gtdbtk/"$bucket"/bacteria.id
done

echo "bacteria names moved"
echo ""

conda activate bactopia 
cd test_gtdbtk
(cd /home/ubuntu/bactopia/datasets/species-specific && ls) > bactopia.data
for i in $(cat bacteria.list)
do
if grep -qw "$i" bactopia.data; then
echo "$i" "already exists" >> bactopia_dataset.log
else
species=$(cat ./taxa/"$i")
bactopia datasets --species "$species" --outdir /home/ubuntu/bactopia/datasets --cpus 16
fi
done
cd ../