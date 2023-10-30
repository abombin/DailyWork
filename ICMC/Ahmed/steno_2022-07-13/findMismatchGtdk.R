setwd("./")

#merged_id<-read.delim("./gtdbtk/output/gtdbtk.bac120.summary.tsv", sep="\t")

merged_id<-readr::read_delim("gtdbtk/output/gtdbtk.bac120.summary.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

mergedIdSep <- data.frame(do.call('rbind', strsplit(as.character(merged_id$classification),';',fixed=TRUE)))
altTax <- data.frame(do.call('rbind', strsplit(as.character(merged_id$`other_related_references(genome_id,species_name,radius,ANI,AF)`),',',fixed=TRUE)))

mergedId<-data.frame(cbind(merged_id$user_genome, mergedIdSep$X7))
mergedAltTax<-data.frame(cbind(merged_id$user_genome, altTax$X2))
colnames(mergedAltTax)<-c("uuid", "AltTax")
mergedAltTax$AltTax<-gsub("s__", "", mergedAltTax$AltTax)

colnames(mergedId)<-c("uuid", "Experiment_Taxa")

mergedId$Experiment_Taxa<-gsub("s__", "", mergedId$Experiment_Taxa)
###
prevID<-readr::read_csv("./metaData/metaDat.csv")
#metaDat <- readr::read_csv("~/ICMC/run_03-14-22/metaData/metaDat.csv")

control1<-prevID[, c(1,5)]
colnames(control1)<-c("uuid", "Control_Taxa")

### part pf the code to be changed after edditing the metadata code
#mergedId$uuid<-paste0("0",mergedId$uuid)
###

control2<-plyr::join(mergedId, control1, by="uuid", type="left", match="first")

control2$Control_Taxa[control2$Control_Taxa=="Clostridium difficile"]<-"Clostridioides difficile"


mismatch <- control2[!(control2$Control_Taxa == control2$Experiment_Taxa), ]



## make list of bacteria that matched

matched<-control2[(control2$Control_Taxa == control2$Experiment_Taxa), ]

matchedSel<-matched[, c(1,3)]
## check if samples mismatched due to misspelling
mismatch$Score<-RecordLinkage::levenshteinSim(mismatch$Control_Taxa, mismatch$Experiment_Taxa)
mismatch$SpellTax <- ifelse((mismatch$Score>0.79), mismatch$Control_Taxa, "NotDecided")
rematched_new<-mismatch[!(mismatch$SpellTax=="NotDecided"),]
rematchedSel<-rematched_new[, c(1,3)]
## check if other closely related genomes would match
notDecided<-mismatch[(mismatch$SpellTax=="NotDecided"),]
altMerged<-plyr::join(notDecided, mergedAltTax,  by="uuid", type="left", match="first")
altMerged$altScore<-RecordLinkage::levenshteinSim(altMerged$Control_Taxa, altMerged$AltTax)
#altMerged$altScore<-stringdist::stringdist(altMerged$Control_Taxa, altMerged$AltTax, method="lv")
altMerged$FinalTax <- ifelse((altMerged$altScore>0.75), altMerged$Control_Taxa, "NotDecided")
altMergedFinal<-altMerged[!(altMerged$FinalTax=="NotDecided"),]
altMergedSel<-altMergedFinal[, c(1,3)]

## bind 
matchedNew<-rbind(matchedSel, rematchedSel, altMergedSel)

setwd("./test_gtdbtk")
## split by taxa
matchedSplit<-split(matchedNew, matchedNew$Control_Taxa)
splitNames<-names(matchedSplit)


for (taxa in splitNames){
  df<- matchedSplit[[taxa]]
  bucket<-tolower(gsub(" ", "-", taxa))
  samples<-data.frame(df$uuid)
  outdir_samples<-paste0("./buckets/",bucket)
  outdir_taxa<-paste0("./taxa/",bucket)
  write.table(samples, file= outdir_samples, col.names = F, row.names = F, quote = F)
  write.table(taxa, file= outdir_taxa, col.names = F, row.names = F, quote = F)
  
}

## working with unmatched samples
notMatchedFinal<-altMerged[(altMerged$FinalTax=="NotDecided"),]
notMatchedSel<-notMatchedFinal[, c(1,4)]
notMatchedSplit<-split(notMatchedSel, notMatchedSel$Experiment_Taxa)
notMatchedSplitNames<-names(notMatchedSplit)

if (length(notMatchedSplit)>0){
  for (i in notMatchedSplitNames){
    df_notMatched<- notMatchedSplit[[i]]
    bucket_notMatched<-tolower(gsub(" ", "-", i))
    samples_notMatched<-data.frame(df_notMatched$uuid)
    outdir_samples_notMatched<-paste0("./buckets/",bucket_notMatched)
    outdir_i_notMatched<-paste0("./i/",bucket_notMatched)
    write.table(samples_notMatched, file= outdir_samples_notMatched, col.names = F, row.names = F, quote = F)
    write.table(i, file= outdir_i_notMatched, col.names = F, row.names = F, quote = F)
  }
} else{
    print("All taxa matched")
}




