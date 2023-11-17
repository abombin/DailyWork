from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import os
from sklearn.metrics import confusion_matrix, roc_auc_score

# redo the predictions with a filter for coverage

curPath = 'IntraSnv_results/ampseq_ConsTest_freq_1_0.csv'
predPath = 'IntraSnv_results/ampseq_comb_derep.csv'
outPath = 'IntraSnv_results/ampseq_ConsTest_freq_1_0_Predictions.csv'

#curPath = 'IntraSnv_results/metaseq_ConsTest_freq_1_0.csv'
#predPath = 'IntraSnv_results/metaseq_comb_derep.csv'
#outPath = 'IntraSnv_results/metaseq_ConsTest_freq_1_0_Predictions.csv'

os.chdir("/home/ubuntu/extraVol/ARVAR/iSNVs")

def getFiltDf(curPath):
    df = pd.read_csv(curPath)
    df = df[df['coverage'] >= 97]
    dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)
    return(dfFilt)

dfFilt = getFiltDf(curPath)
predDf = getFiltDf(predPath)

print(dfFilt['coverage'].min())
print(predDf['coverage'].min())

mask = predDf["Samp_Pos_Ref_Alt"].isin(dfFilt["Samp_Pos_Ref_Alt"])
predDfFilt = predDf[~mask].reset_index(drop = True)

colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq']
colOpt5 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq', 'Sample', 'Sample_AlignPos_Ref_Var']
X_train_1 = dfFilt[colOpt5]
y_train = dfFilt['ConsTest']
#X_train_1, X_test_1, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
X_train = X_train_1[colOpt1]
X_test = predDfFilt[colOpt1]

def getBalancedTrain(X_train, y_train):
    # Combine the features and labels for subsetting
    train_data = pd.concat([X_train, pd.Series(y_train, name="Target")], axis=1)
    # Separate the dataset into class 0 and class 1
    class_0 = train_data[train_data["Target"] == 0]
    class_1 = train_data[train_data["Target"] == 1]
    # Sample an equal number of samples from both classes
    n_samples = min(len(class_0), len(class_1))
    print(n_samples)
    n_samples_0 = int(n_samples * 1.5)
    #print(n_samples_0)
    class_0_subset = class_0.sample(n=n_samples, random_state=42)
    class_1_subset = class_1.sample(n=n_samples, random_state=42)
    # Concatenate the subsets to get a balanced training dataset
    balanced_train_data = pd.concat([class_0_subset, class_1_subset])
    # Split the balanced dataset back into features and labels
    X_train_balanced = balanced_train_data.drop("Target", axis=1)
    y_train_balanced = balanced_train_data["Target"]
    return  X_train_balanced, y_train_balanced

def SumPredict(df, minFreq, maxFreq):
  corPredict = []
  wrongPredict = []
  dfFilt = df[(df["ALLELE.FREQUENCY"] >= minFreq) & (df["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  for i in range(len(dfFilt.index)):
    if dfFilt.loc[i, "Predict_val"] == dfFilt.loc[i, "ConsTest"]:
      corPredict.append(dfFilt.loc[i, "Sample_AlignPos_Ref_Var"])
    else:
      wrongPredict.append(dfFilt.loc[i, "Sample_AlignPos_Ref_Var"])
  numbCorPred = len(corPredict)
  numbWrongPred = len(wrongPredict)
  return numbCorPred, numbWrongPred

def SumPredictPerSamp(df, minFreq, maxFreq):
  snvs = []
  truePos = []
  falsePos = []
  Samples = []
  totalPos = []
  corIdent = []
  wrongIdnet = []
  falseNeg = []
  dfFilt = df[(df["ALLELE.FREQUENCY"] >= minFreq) & (df["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  samples = list(pd.unique(df["Sample"]))
  for sample in samples:
    corPredict = []
    wrongPredict = []
    curTruePos = []
    curFalsePos = []
    curFalseNeg = []
    dfSub = dfFilt[dfFilt["Sample"] == sample].reset_index().drop(["index"], axis =1)
    for i in range(len(dfSub.index)):
      if dfSub.loc[i, "Predict_val"] == dfSub.loc[i, "ConsTest"]:
        corPredict.append(dfSub.loc[i, "Sample_AlignPos_Ref_Var"])
      else:
        wrongPredict.append(dfSub.loc[i, "Sample_AlignPos_Ref_Var"])
        if dfSub.loc[i, "Predict_val"] == 1:
          curFalsePos.append(dfSub.loc[i, "Sample_AlignPos_Ref_Var"])
        elif dfSub.loc[i, "Predict_val"] == 0:
          curFalseNeg.append(dfSub.loc[i, "Sample_AlignPos_Ref_Var"])
    numbCorPred = len(corPredict)
    numbWrongPred = len(wrongPredict)
    dfPos = dfSub[dfSub["ConsTest"] ==1].reset_index().drop(["index"], axis =1)
    curPos = dfPos['Sample_AlignPos_Ref_Var'].nunique()
    curSnvs = dfSub['Sample_AlignPos_Ref_Var'].nunique()
    for i in range(len(dfPos.index)):
      if dfPos.loc[i, "Predict_val"] == dfPos.loc[i, "ConsTest"]:
        curTruePos.append(dfPos.loc[i, "Sample_AlignPos_Ref_Var"])
        
    numbTruePos = len(curTruePos)
    numbFalsePos = len(curFalsePos)
    numbFalseNeg = len(curFalseNeg)
    
    #append values
    snvs.append(curSnvs)
    corIdent.append(numbCorPred)
    wrongIdnet.append(numbWrongPred)
    Samples.append(sample)
    totalPos.append(curPos)
    truePos.append(numbTruePos)
    falsePos.append(numbFalsePos)
    falseNeg.append(numbFalseNeg)
  combDat = pd.DataFrame({"Sample": Samples, "Total_SNVs": snvs, "Correctly_identified": corIdent, "Incorrectly_identified": wrongIdnet , 'False_Positive': falsePos, "False_Negative":falseNeg, 
                          "Total_Positive_Truth": totalPos, "True_Positive_Ident": truePos})
  return combDat

X_train_balanced, y_train_balanced = getBalancedTrain(X_train, y_train)

# Create and fit the Random Forest model on the balanced training dataset
model = RandomForestClassifier(
    criterion='gini',
    n_estimators=3000,
    min_samples_split=2,
    min_samples_leaf=1,
    max_features=None,
    oob_score=True,
    random_state=42,
    n_jobs=-1
)

# redo classification with  balanced complete dataset
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
predDfFilt["ConsTest"] = y_pred

predDfCons = predDfFilt[predDfFilt["ConsTest"] == 1]
predDfCons["Origin"] = "Predict"
overlapCons = dfFilt[dfFilt["ConsTest"] ==1]
overlapCons["Origin"] = "Overlap"

combDat = pd.concat([overlapCons, predDfCons], axis=0, ignore_index=True).reset_index(drop=True)
combDat.head()

combDat.to_csv(outPath, index=False)