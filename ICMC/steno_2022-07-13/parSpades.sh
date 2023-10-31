for i in $(cat ./sample.name)
do 
/home/ubuntu/spades/SPAdes-3.15.4-Linux/bin/spades.py -1 "$i"*val_1.fq.gz -2 "$i"*val_2.fq.gz \
-o spades_out -t 8 --isolate

cp ./spades_out/contigs.fasta ../../gtdbtk/input/"$i".fasta
done

