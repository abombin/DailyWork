/home/ubuntu/nextflow run nf-core/viralrecon -r 2.6.0 \
    # number of cpus to use
    --max_cpus 94 \
    # maximum memory to use
    --max_memory '370.GB' \
    # input file
    --input input.csv \
    # output directory name
    --outdir output_trim_offset \
    --platform illumina \
    # sequencing protocol amplicon or metaseq
    --protocol amplicon \
    # viral genome reference
    --genome 'MN908947.3' \
    # path to kraken database with human reads, if needed
    --kraken2_db ../../references/kraken2-human-db \
    # skip filtering out human reads
    --skip_kraken2 \
    # do not save intemidiate files
    --save_reference false \
    # program to call variants (SNVs)
    --variant_caller 'ivar' \
    # path to the primers file
    --primer_bed ../../references/swift_refv3_primers.bed \
    # primers sufixes, the way for program to know which are the forward and reverse primers, should correspond to the sufixes in the primers .bed file
    --primer_left_suffix "_LEFT" \
    --primer_right_suffix "_RIGHT" \
    # trim additional nucleotides after primer's sequence, currently 0
    --ivar_trim_offset 0 \
    # do not perform de novo genome assembly, only reference based assembly
    --skip_assembly \
    # config file with the versions of nextclade and pangolin 
    -c ../custom.config \
    # profile with dependencies, keep docker as the instance is tuned for it
    -profile docker \
    # name of the docker container
    -with-docker nfcore/virarecon