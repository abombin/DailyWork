cd gtdbtk;


conda activate gtdbtk-2.1.0

mkdir output

gtdbtk classify_wf --genome_dir ./input --out_dir ./output --extension fasta --cpus 10


