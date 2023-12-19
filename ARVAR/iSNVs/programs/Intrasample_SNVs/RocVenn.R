library(ggvenn)
library(ggpubr)
library(MASS)
library(pROC)
library(caret)
library(pscl)
library(car)
library(randomForest)

dir.create("Venn", showWarnings = F)
###
getConsensus <- function(metaSeq, ampSeq, protocol, maxFreq, minFreq, freqCol, filtSteps) {
  if (filtSteps == 1) {
    metaseq <- read.csv(metaSeq)
    ampseq <- read.csv(ampSeq)
    
    metaseq =  metaseq[metaseq$coverage >= 97,]
    ampseq =   ampseq[ampseq$coverage >= 97,]
    
    metaseqFilt <- metaseq[metaseq[[freqCol]] >= minFreq & metaseq[[freqCol]] <= maxFreq, ]
    ampseqFilt <- ampseq[ampseq[[freqCol]] >= minFreq & ampseq[[freqCol]] <= maxFreq, ]
    
    ampseqFilt = ampseqFilt[ampseqFilt$Sample%in%metaseqFilt$Sample,]
    metaseqFilt = metaseqFilt[metaseqFilt$Sample%in%ampseqFilt$Sample,]
    if (protocol == "ampseq") {
      # dataframe with stats
      targDf = ampseqFilt
      # data to check agains
      refDf = metaseqFilt
      targSnv <- unique(refDf$Samp_Pos_Ref_Alt)
    } else if (protocol == "metaseq") {
      # dataframe with stats
      targDf = metaseqFilt
      # data to check agains
      refDf = ampseqFilt
      targSnv <- unique(refDf$Samp_Pos_Ref_Alt)
    }
    ConsTest <- numeric()
    for (i in 1:nrow(targDf)) {
      curSnv =  targDf[i, "Samp_Pos_Ref_Alt"]
      if (curSnv%in%targSnv){
        curCons = 1
      } else {
        curCons = 0
      }
      ConsTest = c(ConsTest, curCons)
    }
    targDf$ConsTest <- ConsTest
    return(targDf)
  } else if (filtSteps == 2 & protocol == "metaseq") {
    df = read.csv(metaSeq) 
    dfFilt = df[df[[freqCol]] >= minFreq & df[[freqCol]] <= maxFreq, ]
    return(dfFilt)
  } else if (filtSteps == 2 & protocol == "ampseq") {
    df = read.csv(ampSeq) 
    dfFilt = df[df[[freqCol]] >= minFreq & df[[freqCol]] <= maxFreq, ]
    return(dfFilt)
  } 
}


runRoc = function(df, protocol, freqCol, splitPerc) {
  dfFilt = df[!is.na(df$Var_Al_RelPos),]
  
  set.seed(42)
  train_idx <- createDataPartition(dfFilt$ConsTest, p = splitPerc, list = FALSE)
  train_data <- dfFilt[train_idx, ]
  test_data <- dfFilt[-train_idx, ]
  
  if (protocol == "metaseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + I(meandepth^2)" ))
    aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    
    # glm_formula <- as.formula(paste("ConsTest ~", freqCol, "+ STRAND.BIAS + QUAL + Var_Al_RelPos + meandepth" ))
    # aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 1000)
    
  } else if (protocol == "ampseq") {
    
    glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + I(meandepth^2)"))
    aucModel <- glm(formula = glm_formula, data = train_data, family = "binomial")
    
    # glm_formula <- as.formula(paste("ConsTest ~", freqCol, " + QUAL + Var_Al_RelPos  + meandepth"))
    # aucModel <- randomForest(formula = glm_formula, data = train_data, ntree = 1000)
    
  }
  
  probs <- predict(aucModel, newdata = test_data, type = "response")
  roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)
  
  # actual predictions and testing
  threshold <- 0.5
  test_data$PredictedConsTest <- ifelse(probs >= threshold, 1, 0)
  equal_count <- sum(test_data$ConsTest == test_data$PredictedConsTest)
  # Count the number of times values in Column1 are NOT equal to Column2
  not_equal_count <- sum(test_data$ConsTest != test_data$PredictedConsTest)
  sumCorrect = equal_count / (equal_count + not_equal_count) * 100
  print(roc_obj$auc)
  print(sumCorrect )
}


getVenn<-function(metaseq, ampseq, filtSteps, maxFreq, minFreq){
  
  ampSnps<-ampseq$Samp_Pos_Ref_Alt
  # process meta
  
  metaSnps<-metaseq$Samp_Pos_Ref_Alt
  print(length(metaSnps))
  print(length(ampSnps))
  print(length(metaSnps[metaSnps%in%ampSnps]))
  print(length(metaSnps[metaSnps%in%ampSnps]) / (length(metaSnps) + length(ampSnps)) * 100 )
  # make a list
  snps<-list('Amplicon'=ampSnps, 'Metagenomic'=metaSnps)
  curTitle = paste0("FiltSteps_", filtSteps, "_MaxFreq_", maxFreq, "_minFreq_", minFreq)
  # plot
  venPlot<-ggvenn(snps, text_size=5, set_name_size=8)+
    ggtitle(curTitle)+
    theme(text = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5, size = 24))
  #ggsave(filename =  paste0("Venn/Intra_", curTitle,'.jpeg'), plot = venPlot, width = 9, height = 9, units = 'in', dpi = 600, device = 'jpeg')
}


# no frequency filtering
metaseq = getConsensus(metaSeq="IntraSnv_results/metaseq_comb_derep_overlap.csv", ampSeq="IntraSnv_results/ampseq_comb_derep_overlap.csv", 
                       protocol="metaseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY", filtSteps=1)
ampseq = getConsensus(metaSeq="IntraSnv_results/metaseq_comb_derep_overlap.csv", ampSeq="IntraSnv_results/ampseq_comb_derep_overlap.csv", 
                      protocol="ampseq", maxFreq=1, minFreq=0, freqCol="ALLELE.FREQUENCY", filtSteps=1)



getVenn(metaseq=metaseq, ampseq=ampseq, filtSteps=1, maxFreq = 1, minFreq = 0)


runRoc(df=metaseq, protocol="metaseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)
runRoc(df=ampseq, protocol="ampseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)


## get the coverage comparisons
ampStats = unique(ampseq[, c("Sample", "coverage", "endpos", "OrigName")])
ampStats = ampStats[order(ampStats$Sample),]
metaStats = unique(metaseq[, c("Sample", "coverage", "endpos", "OrigName")])
metaStats = metaStats[order(metaStats$Sample),]
identical(ampStats$Sample, metaStats$Sample)
identical(ampStats$coverage, metaStats$coverage)
colnames(metaStats)[2:4] = c("coverage_meta", "endpos_meta", "OrigName_meta")
combStats = plyr::join(ampStats, metaStats, by = "Sample", type = "left", match = "all")
identical(combStats$endpos, combStats$endpos_meta)

write.csv(combStats, "ampseq_metaseq_overlap_combStats_2.csv", row.names = F)

# try without indels
metaseqSub = metaseq[nchar(metaseq$REF.NT) == 1 & nchar(metaseq$VAR.NT) == 1, ]
ampSub = ampseq[nchar(ampseq$REF.NT) == 1 & nchar(ampseq$VAR.NT) == 1, ]

runRoc(df=metaseqSub, protocol="metaseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)
runRoc(df=ampSub, protocol="ampseq", freqCol="ALLELE.FREQUENCY", splitPerc=0.7)

getVenn(metaseq=metaseqSub, ampseq=ampSub, filtSteps=1, maxFreq = 1, minFreq = 0)

# make histograms for frequency distributions
# make tables SNVs per Sample at freq 0, 0.01 and 0.02