conda activate bactopia-dev

bactopia --samples bactopiaInput.tsv --species "Stenotrophomonas maltophilia" --outdir bactopia_output --max_cpus 8

rm -rf ./panOut
rm -rf work
bactopia --wf pangenome --bactopia ./bactopia_output --outdir ./panOut --max_cpus 16 --skip_phylogeny --skip_recombination --run_name 'all_samples'

rm -rf work

