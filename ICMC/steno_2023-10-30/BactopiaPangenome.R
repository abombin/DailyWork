system("bactopia --samples bactopiaInput.tsv --species 'Stenotrophomonas maltophilia' --outdir bactopia_output --max_cpus 8")

system("rm -rf ./panOut")
system("rm -rf work")
system("bactopia --wf pangenome --bactopia ./bactopia_output --outdir ./panOut --max_cpus 16 --skip_phylogeny --skip_recombination --run_name 'all_samples'")

system("rm -rf work")

curDirs = list.files(paste0("panOut/bactopia-runs"), full.names = T)
# change next line when run with normal bactopia runs
for (curDir in curDirs) {
  setwd(curDir)
  iqCommand<-paste0('iqtree -s core-genome.aln.gz -m TEST -nt 30 -bb 1000')
  system(iqCommand)
  setwd("../../../")
}