cd processPar

ls -d * | parallel ' cd {} && /home/ubuntu/spades/SPAdes-3.15.5-Linux/bin/spades.py -1 sample_val_1.fq.gz -2 sample_val_2.fq.gz -o spades_out -t 8 --isolate'