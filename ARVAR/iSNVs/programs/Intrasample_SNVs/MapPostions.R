library(foreach)
library(doParallel)

combStats = read.csv("ampseq_metaseq_overlap_combStats.csv")

dir.create("Overlap_Pos_mapping/Comb_fasta/", recursive = T, showWarnings = F)

makeAlignment = function(combStats, i) {
  ampSample = combStats$OrigName[i]
  metaSample = combStats$OrigName_meta[i]
  sampName = combStats$Sample[i]
  catFile = paste0("Overlap_Pos_mapping/Comb_fasta/",   sampName, ".fasta")
  ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
  metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
  cmd_str = paste("cat references/MN908947.3.fna", ampPath, metaPath, ">", catFile, sep = " ")
  system(cmd_str)
  cmd_mafft = paste0("mafft --maxiterate 1000 --thread 4 --globalpair ",  catFile, " > Overlap_Pos_mapping/", sampName, "_align.fasta")
  system(cmd_mafft)
}


dir.create("Overlap_Pos_mapping/Indecies/", recursive = T, showWarnings = F)
dir.create("Overlap_Pos_mapping/Indecies/Ampseq", recursive = T, showWarnings = F)
dir.create("Overlap_Pos_mapping/Indecies/Metaseq", recursive = T, showWarnings = F)


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

compileIndex = function(combStats, i, makeCheck) {
  ampSample = combStats$OrigName[i]
  metaSample = combStats$OrigName_meta[i]
  sampName = combStats$Sample[i]
  ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
  metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
  alignmentPath = paste0("Overlap_Pos_mapping/", sampName, "_align.fasta")
  
  # read fasta
  curAlignment = phylotools::read.fasta(alignmentPath)
  
  # split strings to vectors
  ampAlignVect = toupper(strsplit(curAlignment[2,2], "")[[1]])
  metaAlignVect = toupper(strsplit(curAlignment[3,2], "")[[1]])
  
  ampIndex = makeIndex(ampAlignVect)
  metaIndex = makeIndex(metaAlignVect)
  
  ampIndex$Sample = sampName
  ampIndex$OrigName =  ampSample 
  metaIndex$Sample = sampName
  metaIndex$OrigName =  metaSample 
  
  # next step is to write each index table in a separate folder for metaseq and ampseq, filename based on the original name
  # in SNV tables add Alignment position and make new SNP as Sample_AlignPos_VarAllele
  outAmp = paste0("Overlap_Pos_mapping/Indecies/Ampseq/", ampSample, ".csv")
  outMeta = paste0("Overlap_Pos_mapping/Indecies/Metaseq/", metaSample, ".csv")
  
  write.csv(ampIndex, outAmp, row.names = F)
  write.csv(metaIndex, outMeta, row.names = F)
  
  if (makeCheck==T) {
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

#makeAlignment(combStats=combStats, i=57)
#compileIndex(combStats=combStats, i=1, makeCheck=T)

cl = makeCluster(7, type = "FORK")
registerDoParallel(cl)
combAlig = foreach(i=1:nrow(combStats))%dopar% {
  makeAlignment(combStats=combStats, i=i)
}
stopCluster(cl)


cl = makeCluster(28, type = "FORK")
registerDoParallel(cl)
combCompile = foreach(i=1:nrow(combStats))%dopar% {
  compileIndex(combStats=combStats, i=i, makeCheck=T)
}
stopCluster(cl)

combCompChar = unlist(combCompile)

write.table(combCompChar, "Overlap_index.log", row.names = F, col.names = F, quote = F)

