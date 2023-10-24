#!/usr/bin/env bash

cat references/MN908947.3.fna IntraSnv_ampseq_overlap/EHC-C19-3440H_S44_L001/reference.fa > comb_amp_EHC-C19-3440H_S44_L001.fasta

mafft --maxiterate 1000 --thread 4 --globalpair comb_amp_EHC-C19-3440H_S44_L001.fasta > amp_alignment_EHC-C19-3440H_S44_L001.fasta

blastn -query IntraSnv_ampseq_overlap/EHC-C19-3440H_S44_L001/reference.fa -subject references/MN908947.3.fna > amp_alignment_EHC-C19-3440H_S44_L001_blast.fasta

mafft --addfragments IntraSnv_ampseq_overlap/EHC-C19-3440H_S44_L001/reference.fa --keeplength references/MN908947.3.fna > aligned_sequences.fasta