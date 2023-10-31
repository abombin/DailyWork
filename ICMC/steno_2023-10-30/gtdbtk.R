# conda activate gtdbtk-2.1.0


runGtdb = function() {
  dir.create('gtdbtk/output', showWarnings = F, recursive = T)
  cmd_str = "gtdbtk classify_wf --genome_dir gtdbtk/input --out_dir gtdbtk/output --extension fasta --cpus 20 --skip_ani_screen"
  system(cmd_str)
}


runGtdb()

gtdb<-readr::read_delim("gtdbtk/output/gtdbtk.bac120.summary.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# separate main gtdbtk classification
processGtdb<-function(x){
  sepNames<-data.frame(do.call('rbind', strsplit(as.character(x$classification),';',fixed=TRUE)))
  gtdbName<-data.frame(cbind(x$user_genome, sepNames$X7, sepNames$X6))
  colnames(gtdbName)<-c("uuid", "Experiment_taxa", "Experiment_genus")
  gtdbName$Experiment_taxa<-gsub("s__", "", gtdbName$Experiment_taxa)
  gtdbName$Experiment_taxa[gtdbName$Experiment_taxa==""]<-'unknown'
  return(gtdbName)
}

# get main species from gtdb classification
gtdb_result<-processGtdb(gtdb)


ncbiDict = read.delim("/home/ubuntu/gtdbtk/conversion/dictionary/gtdbtkNcbiDictSpecies.tsv")

unique(ncbiDict[grepl("Stenotrophomonas sp003028475", ncbiDict$gtdb_species),])