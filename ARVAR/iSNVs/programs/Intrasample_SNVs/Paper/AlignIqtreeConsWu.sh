cat references/MN908947.3.fna > Paper/allAlignment/Comb_fasta/allCons.fasta

echo >> Paper/allAlignment/Comb_fasta/allCons.fasta

cat Consens_fasta/Ampseq/*.fa >> Paper/allAlignment/Comb_fasta/allCons.fasta

echo >> Paper/allAlignment/Comb_fasta/allCons.fasta

cat Consens_fasta/Metaseq/*.fa >> Paper/allAlignment/Comb_fasta/allCons.fasta

mafft --auto --thread 32 Paper/allAlignment/Comb_fasta/allCons.fasta > Paper/allAlignment/allCons.align

cd /home/ubuntu/extraVol/ARVAR/iSNVs/Paper/allAlignment/

# added row names with __ separator to alignment names  
time iqtree -s allCons_edit.align -m TEST -nt 30 -bb 1000

