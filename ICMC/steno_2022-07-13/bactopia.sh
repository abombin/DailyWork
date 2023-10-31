conda activate bactopia

(cd ./bactopia_gtdbtk && ls) > bacteria.list
cd bactopia_gtdbtk
for i in $(cat ../bacteria.list); do cd ./"$i"
mkdir output
organism=$(cat bacteria.id)

bactopia prepare ./input > fastqs.txt

bactopia --fastqs fastqs.txt \
	--datasets /home/ubuntu/bactopia/datasets \
	--species "$organism" \
	--coverage 100 \
	--genome_size median \
	--outdir ./output/ \
	--max_cpus 8

cd ../
done
