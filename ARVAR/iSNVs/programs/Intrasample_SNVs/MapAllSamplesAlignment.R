library(foreach)
library(doParallel)
library(phylotools)

combStats = read.csv("ampseq_metaseq_overlap_combStats.csv")

dir.create("Overlap_Pos_mapping/Comb_fasta/", recursive = T, showWarnings = F)

makeAlignment = function(combStats) {
  #system("cat references/MN908947.3.fna > Overlap_Pos_mapping/Comb_fasta/allSequences.fasta")
  seqType = character()
  seqNames = character()
  for ( i in 1:nrow(combStats) ) {
    ampSample = combStats$OrigName[i]
    metaSample = combStats$OrigName_meta[i]
    sampName = combStats$Sample[i]
    seqNames= c(seqNames, ampSample, metaSample)
    seqType = c(seqType, "ampseq", "metaseq")
    
    catFile = "Overlap_Pos_mapping/Comb_fasta/allSequences.fasta"
    ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
    metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
    
    cmd_str = paste("cat ", ampPath, metaPath, ">>", catFile, sep = " ")
    #system(cmd_str)
  }
  cmd_mafft = paste0("mafft --maxiterate 1000 --thread 29 --auto ",  catFile, " > Overlap_Pos_mapping/All_samples_align.fasta")
  system(cmd_mafft)
  
  combDat = data.frame(seqNames, seqType)
  write.csv(combDat, "All_overlap_samples_alignment_order.csv", row.names = F)
}

# run alignment function
makeAlignment(combStats)

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

alignOrder = read.csv("All_overlap_samples_alignment_order.csv")
alignFasta = read.fasta("Overlap_Pos_mapping/All_samples_align.fasta")
alignFasta = alignFasta[2:nrow(alignFasta),]
alignFasta$seq.name = gsub(" MN908947.3", "", alignFasta$seq.name)
identical(alignFasta$seq.name, alignOrder$seqNames)
alignFasta$SeqType = alignOrder$seqType

compileIndex = function(i, makeCheck, alignFasta) {
  print(i)
  ampAlignment = alignFasta[alignFasta$SeqType == "ampseq",]
  rownames(ampAlignment) = NULL
  metaAlignment=alignFasta[alignFasta$SeqType == "metaseq",]
  rownames( metaAlignment) = NULL
  
  # split strings to vectors
  ampAlignVect = toupper(strsplit(ampAlignment[i,2], "")[[1]])
  metaAlignVect = toupper(strsplit(metaAlignment[i,2], "")[[1]])
  
  ampIndex = makeIndex(ampAlignVect)
  metaIndex = makeIndex(metaAlignVect)
  
  ampIndex$OrigName = ampAlignment[i,1]
  metaIndex$OrigName = metaAlignment[i,1]
  
  ampSample = ampAlignment[i,1]
  metaSample = metaAlignment[i,1]
  
  outAmp = paste0("Overlap_Pos_mapping/AllSampAlign_Indecies/Ampseq/", ampSample, ".csv")
  outMeta = paste0("Overlap_Pos_mapping/AllSampAlign_Indecies/Metaseq/", metaSample, ".csv")
  
  write.csv(ampIndex, outAmp, row.names = F)
  write.csv(metaIndex, outMeta, row.names = F)
  
  if (makeCheck==T) {
    ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
    metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
    
    ampFasta = phylotools::read.fasta(ampPath)
    metaFasta = phylotools::read.fasta(metaPath)
    
    ampVector = toupper(strsplit(ampFasta[1,2], "")[[1]])
    metaVector = toupper(strsplit(metaFasta[1,2], "")[[1]])
    
    checkAmp = checkIndex(combIndex=ampIndex, targVec=ampVector, alignTargVec=ampAlignVect)
    checkMeta = checkIndex(combIndex=metaIndex, targVec=metaVector, alignTargVec=metaAlignVect)
    
    print(paste0("Ampseq sample ", ampSample, " failed rows ", length(checkAmp)))
    print(paste0("Metaseq sample ", metaSample, " failed rows ", length(checkMeta)))
    
    checkLog = c(paste0("Ampseq sample ", ampSample, " failed rows ", length(checkAmp)), paste0("Metaseq sample ", metaSample, " failed rows ", length(checkMeta)))
    return(checkLog)
  }

}

dir.create("Overlap_Pos_mapping/AllSampAlign_Indecies/Ampseq/", recursive = T, showWarnings = F)
dir.create("Overlap_Pos_mapping/AllSampAlign_Indecies/Metaseq/",  recursive = T, showWarnings = F)

cl = makeCluster(14, type = "FORK")
registerDoParallel(cl)
combCompile = foreach(i=1:(nrow(alignFasta))/2)%dopar% {
  compileIndex(i=i, makeCheck=T, alignFasta=alignFasta)
}
stopCluster(cl)

combCompChar = unlist(combCompile)

write.table(combCompChar, "Overlap_index_allSampAlign.log", row.names = F, col.names = F, quote = F)

combCompile = character()
for ( i in 1:nrow(alignFasta)) {
  curCompile = compileIndex(i=i, makeCheck=T, alignFasta=alignFasta)
  combCompile = c(combCompile, curCompile)
}
