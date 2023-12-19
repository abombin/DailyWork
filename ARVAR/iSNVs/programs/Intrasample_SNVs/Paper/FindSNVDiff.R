library(phylotools)

fasta = read.fasta("Paper/allAlignment/allCons.align")

fasta = read.fasta("Paper/allAlignment/allSamples.align")

sumFreq = data.frame(table(fasta$seq.name))

fastaDup = fasta[fasta$seq.name == "EHC-C19-1463I_S15_L001 MN908947.3",]

dat2fasta(fastaDup, "EHC-C19-1463I_S15_L001_ConsWu.align")
dat2fasta(fastaDup, "EHC-C19-1463I_S15_L001_Int.align")

system("snp-dists/snp-dists EHC-C19-1463I_S15_L001_ConsWu.align > EHC-C19-1463I_S15_L001_ConsWu.dist")
system("snp-dists/snp-dists EHC-C19-1463I_S15_L001_Int.align > EHC-C19-1463I_S15_L001_Int.dist")


distMat =  read.delim("~/extraVol/ARVAR/iSNVs/EHC-C19-1463I_S15_L001_ConsWu.dist")
distMat = read.delim("~/extraVol/ARVAR/iSNVs/EHC-C19-1463I_S15_L001_Int.dist")

corChar = c("a", "t", "g", "c")


samp1 =  strsplit(fastaDup$seq.text[1], split ="")[[1]]
samp2 = strsplit(fastaDup$seq.text[2], split ="")[[1]]

difBases = character()
for ( i in 1:nchar(fastaDup$seq.text[1])) {
  base1 = samp1[i]
  base2 = samp2[i]
  if (base1%in%corChar &  base2%in%corChar & base1!=base2) {
    print(i)
    print(base1)
    print(base2)
  }
}