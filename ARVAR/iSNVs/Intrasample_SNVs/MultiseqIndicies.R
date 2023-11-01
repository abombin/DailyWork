library(phylotools)
library(foreach)
library(doParallel)

# make index
makeIndex = function(alignTargVec) {
  Alignment_Pos = numeric()
  Original_Pos = numeric()
  alignBases = character()
  curOrigIndex = 0
  for ( i in 1:length(alignTargVec) ) {
    curBase = alignTargVec[i]
    if ( curBase != "-" ) {
      Alignment_Pos = c(Alignment_Pos, i)
      curOrigIndex = curOrigIndex + 1
      Original_Pos = c(Original_Pos, curOrigIndex)
      #alignBases = c(alignBases,  curBase)
    }
  }
  combIndex = data.frame(Original_Pos, Alignment_Pos)
  return(combIndex)
}

checkIndex = function(combIndex, targVec, alignTargVec) {
  combChecks = character()
  for ( i in 1:nrow(combIndex) ) {
    curOriginalPos = targVec[combIndex$Original_Pos[i]]
    curAlignPos = alignTargVec[combIndex$Alignment_Pos[i]]
    if (curOriginalPos !=  curAlignPos ) {
      combChecks = c(combChecks, i)
    }
  }
  return(combChecks)
}

# 


# start 
metaAlign = read.fasta("/home/ubuntu/extraVol/ARVAR/iSNVs/All_Samples_Mapping/metaseq/All_samples_align.fasta")
ampseqAlign = read.fasta("/home/ubuntu/extraVol/ARVAR/iSNVs/All_Samples_Mapping/ampseq/All_samples_align.fasta")

metaFasta = read.fasta('All_Samples_Mapping/metaseq/Comb_fasta/allSequences.fasta')
ampFasta = read.fasta('All_Samples_Mapping/ampseq/Comb_fasta/allSequences.fasta')

alignment = read.fasta("All_Samples_Mapping/All/All_samples_align.fasta")
curAlignment = alignment[!alignment$seq.name == "MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome",]

meta_vector <- rep("metaseq", times = 399)
amp_vector = rep("ampseq", times = 712)
combVect = c(amp_vector, meta_vector)

curAlignment$SeqType = combVect


# ampseqAlign[2:nrow(ampseqAlign), 1]
# curAlignment[1:712, 1]

identical(ampseqAlign[2:nrow(ampseqAlign), 1], curAlignment[1:712, 1])


curAlignment$seq.name = gsub(" MN908947.3", "", curAlignment$seq.name)

dir.create("All_Samples_Mapping/AllSampAlign_Indecies/ampseq/", recursive = T)
dir.create("All_Samples_Mapping/AllSampAlign_Indecies/metaseq/", recursive = T)


compileIndex = function(i, makeCheck, alignFasta, protocol) {
  
  # split strings to vectors
  
  
  AlignVect = toupper(strsplit(alignFasta[i,2], "")[[1]])
  
  AlignIndex = makeIndex(AlignVect)
  
  AlignIndex$OrigName = alignFasta[i,1]
  
  SampleName = alignFasta[i,1]
  
  outFile= paste0("All_Samples_Mapping/AllSampAlign_Indecies/",protocol,"/", SampleName, ".csv")

  write.csv(AlignIndex,  outFile, row.names = F)
  
  if (makeCheck==T) {
    
    if (protocol == "ampseq") {
      ampPath = paste0("IntraSnv_ampseq_all/",  SampleName, "/reference.fa")
      ampFasta = phylotools::read.fasta(ampPath)
      ampVector = toupper(strsplit(ampFasta[1,2], "")[[1]])
      checkSample = checkIndex(combIndex=AlignIndex, targVec=ampVector, alignTargVec=AlignVect)
      
    } else if (protocol == "metaseq") {
      metaPath = paste0("IntraSnv_metaseq_all/",  SampleName, "/reference.fa")
      metaFasta = phylotools::read.fasta(metaPath)
      metaVector = toupper(strsplit(metaFasta[1,2], "")[[1]])
      checkSample = checkIndex(combIndex=AlignIndex, targVec=metaVector, alignTargVec=AlignVect)
    }
    
    
    print(paste0("Sample ", SampleName, " failed rows ", length(checkSample)))
    
    checkLog = c(paste0("Sample ", SampleName, " failed rows ", length(checkSample)))
    return(checkLog)
  }
  
}


ampAlignment = curAlignment[curAlignment$SeqType == "ampseq",]
rownames(ampAlignment) = NULL

metaAlignment=curAlignment[curAlignment$SeqType == "metaseq",]
rownames( metaAlignment) = NULL

cl = makeCluster(15, type = "FORK")
registerDoParallel(cl)
combCompile = foreach(i=1:(nrow(ampAlignment)))%dopar% {
  compileIndex(i=i, makeCheck=T, alignFasta=ampAlignment, protocol = "ampseq")
}
stopCluster(cl)

combCompChar = unlist(combCompile)
write.table(combCompChar, "AllSamplesAlignmentIndex_ampseq.log", row.names = F, col.names = F, quote = F)


cl = makeCluster(15, type = "FORK")
registerDoParallel(cl)
combCompile = foreach(i=1:(nrow(metaAlignment)))%dopar% {
  compileIndex(i=i, makeCheck=T, alignFasta=metaAlignment, protocol = "metaseq")
}
stopCluster(cl)

combCompChar = unlist(combCompile)
write.table(combCompChar, "AllSamplesAlignmentIndex_ampseq.log", row.names = F, col.names = F, quote = F)