
refFasta = phylotools::read.fasta("references/MN908947.3.fna")
refFasta[1,2]

targFasta = phylotools::read.fasta("IntraSnv_ampseq_overlap/EHC-C19-3440H_S44_L001/reference.fa")
targFasta[1,2]
nchar(targFasta[1,2])

curAlignment = phylotools::read.fasta('amp_alignment_EHC-C19-3440H_S44_L001.fasta')
curAlignment[1,2]
nchar(curAlignment[1,2])

curAlignment[2,2]
nchar(curAlignment[2,2])



alignTargVec = toupper(strsplit(curAlignment[2,2], "")[[1]])
targVec = toupper(strsplit(targFasta[1,2], "")[[1]])


alignIndex = numeric()
originalIndex = numeric()
alignBases = character()

curOrigIndex = 0
for ( i in 1:length(alignTargVec) ) {
  curBase = alignTargVec[i]
  
  if ( curBase != "-" ) {
    alignIndex = c(alignIndex, i)
    curOrigIndex = curOrigIndex + 1
    originalIndex = c(originalIndex, curOrigIndex)
    alignBases = c(alignBases,  curBase)
  }
}

combIndex = data.frame(originalIndex, alignIndex)

alignTargVec[6516]
targVec[6513]

length(alignBases)
length(targVec)
nrow(combIndex)

DiffInd = combIndex[combIndex$originalIndex != combIndex$alignIndex,]

colnames(combIndex)

# make a check that basses are correct 

combChecks = character()

for ( i in 1:nrow(combIndex) ) {
  curOriginalPos = targVec[combIndex$originalIndex[i]]
  curAlignPos = alignTargVec[combIndex$alignIndex[i]]
  if (curOriginalPos !=  curAlignPos ) {
    combChecks = c(combChecks, i)
  }
}