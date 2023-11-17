from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import os
from sklearn.metrics import confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

os.chdir("/home/ubuntu/extraVol/ARVAR/iSNVs")

df = pd.read_csv('IntraSnv_results/metaseq_ConsTest_freq_1_0.csv')
#plotname = "Metaseq_auc_plot.png"
#df = pd.read_csv('IntraSnv_results/ampseq_ConsTest_freq_1_0.csv')


#df = pd.read_csv('IntraSnv_ampseq_overlap/ampseq_metaseq_overlap_97_allFreq.csv')

mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)


# make frequency filtering then apply is.in to filter out samples that do not overlap
#colOpt1 = ["ALLELE.FREQUENCY", "STRAND.BIAS" , "DEPTH", "QUAL", "Var_Al_RelPos", "Ref_Al_RelPos", "meandepth", "coverage",  "meanmapq", "meanbaseq"]
#colOpt5 = ["ALLELE.FREQUENCY", "STRAND.BIAS" , "DEPTH", "QUAL", "Var_Al_RelPos", "Ref_Al_RelPos", "meandepth", "coverage",  "meanmapq", "meanbaseq", 'Sample', 'Sample_AlignPos_Ref_Var']
#colOpt1 = ["ALLELE.FREQUENCY", "STRAND.BIAS" , "QUAL", "meandepth", "coverage",  "meanmapq", "meanbaseq"]
#colOpt1 = ['ALLELE.FREQUENCY', 'DEPTH', 'Var_Al_RelPos', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq', 'STRAND.BIAS', 'QUAL']


colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq']
colOpt5 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq', 'Sample', 'Sample_AlignPos_Ref_Var']

X = dfFilt[colOpt5]
y = dfFilt['ConsTest']

varNa = X[X.isna().any(axis=1)]
respNa = y[y.isna()]

X_train_1, X_test_1, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

X_train = X_train_1[colOpt1]
X_test = X_test_1[colOpt1]

model = RandomForestClassifier(criterion='gini', 
                             n_estimators=3000,
                             min_samples_split=2,
                             min_samples_leaf=1,
                             max_features=None,
                             oob_score=True,
                             random_state=42,
                             n_jobs=-1)
model.fit(X_train, y_train)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")

print("%.4f" % model.oob_score_)
# show variables importance 
pd.concat((pd.DataFrame(X_train.columns, columns = ['variable']), 
           pd.DataFrame(model.feature_importances_, columns = ['importance'])), 
          axis = 1).sort_values(by='importance', ascending = False)[:20]


varImportance = pd.concat((pd.DataFrame(X_train.columns, columns = ['variable']), 
           pd.DataFrame(model.feature_importances_, columns = ['importance'])), 
          axis = 1).sort_values(by='importance', ascending = False)[:20]

varImportance["AUC"] = auc_score

#varImportance.to_csv("IntraSnv_results/SumModels/Ampseq_VarImportance.csv", index=False)
#varImportance.to_csv("IntraSnv_results/SumModels/Metaseq_VarImportance.csv", index=False)


y_pred = model.predict(X_test)
cm = confusion_matrix(y_test, y_pred)

# Calculate percentages
print(cm)

# Extract the true negative, false positive, false negative, and true positive values
tn, fp, fn, tp = cm.ravel()

# Calculate class error rates
class_0_error = fp / (tn + fp)  # Class 0 (ConsTest = 0)
class_1_error = fn / (fn + tp)  # Class 1 (ConsTest = 1)

# Print the class error rates
print(f"Class 0 Error Rate: {class_0_error:.4f}")
print(f"Class 1 Error Rate: {class_1_error:.4f}")




# make data with Ids and survival
X_test_1["Predict_val"] = y_pred
X_test_1['ConsTest'] = y_test

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

numbCorPred, numbWrongPred = SumPredict(df = X_test_1, minFreq = 0, maxFreq = 1)

numbCorPred / (numbCorPred + numbWrongPred) * 100


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

sumTest = SumPredictPerSamp(df = X_test_1, minFreq = 0, maxFreq = 1)

sumTest["True_Positive_Ident"].sum() / sumTest["Total_Positive_Truth"].sum() * 100

#sumTest.to_csv("IntraSnv_results/ampseq_ConsTest_freq_1_0_predictSum.csv", index=False)

print(sumTest["True_Positive_Ident"].sum() / sumTest["Total_Positive_Truth"].sum() * 100)

print(sumTest['Correctly_identified'].sum() / sumTest['Total_SNVs'].sum() * 100)

print(f"AUC Score: {auc_score}")

# make plots 

def plot_auc_curve(y_test, y_pred_proba, plotname):
    # Create a new figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))  # Adjust width and height as desired

    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)

    # Clear the previous plot
    ax.clear()

    # Plotting the AUC curve
    ax.plot(fpr, tpr, label='AUC Curve (AUC = %0.2f)' % auc_score)
    ax.plot([0, 1], [0, 1], 'k--', label='Random')
    ax.set_xlabel('False Positive Rate (FPR)', fontsize=18)
    ax.set_ylabel('True Positive Rate (TPR)', fontsize=18)
    ax.set_title('Receiver Operating Characteristic (ROC)', fontsize=18)
    ax.legend(loc='lower right',  fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)

    # Save the plot with adjusted width, height, and DPI
    plt.savefig(plotname, dpi=300, bbox_inches='tight')

    # Display the plot
    plt.show()
    
plot_auc_curve(y_test = y_test, y_pred_proba = y_pred_proba, plotname = plotname)
