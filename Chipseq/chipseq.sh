/home/flyhunter/nextflow run nf-core/chipseq -r 2.0.0 \
--input input.csv --outdir output_bowtie2 \
--fasta /home/flyhunter/Wang/data/reference/refdata-gex-mm10-2020-A/fasta/genome.fa \
--gtf /home/flyhunter/Wang/data/reference/refdata-gex-mm10-2020-A/genes/genes.gtf \
--blacklist /home/flyhunter/Kai/Chipseq/references/mm10-blacklist.v2.bed \
--narrow_peak \
--aligner 'bowtie2' \
--macs_fdr 0.05 \
--macs_gsize 1870000000 \
--max_cpus 14 \
--max_memory '58.GB' \
--save_reference \
-profile docker \
-with-docker nfcore/chipseq