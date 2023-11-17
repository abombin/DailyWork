from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd



minFreq = 0.02
maxFreq = 1

def getConsensus(metaSeq, ampSeq, minFreq, maxFreq):
  metaseq = pd.read_csv(metaSeq)
  ampseq = pd.read_csv(ampSeq)
  metaseqFilt = metaseq[(metaseq["ALLELE.FREQUENCY"] >= minFreq) & (metaseq["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  ampseqFilt =  ampseq[(ampseq["ALLELE.FREQUENCY"] >= minFreq) & (ampseq["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  targSnv = metaseqFilt["Samp_Pos_Ref_Alt"].to_list()
  ConsTest = []
  for i in range(len(ampseqFilt.index)):
    curSnv = ampseqFilt.loc[i, "Samp_Pos_Ref_Alt"]
    if curSnv in targSnv:
      curCons = 1
    else:
      curCons = 0
    ConsTest.append(curCons)
  ampseqFilt["ConsTest"] = ConsTest
  return(ampseqFilt)

# increasing lower frequency increase accuracy 
# colOption1, minFreq at least 0.02 and ampSeqIvar have the bset performance   
dfFilt = getConsensus(metaSeq = "snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_97_v2.csv", ampSeq = "snvs_comb_res/ampseq_overlap_comb_derep_covFilt_97_v2.csv", minFreq = minFreq, maxFreq = maxFreq)
dfRef = dfFilt[dfFilt["ConsTest"] ==1].reset_index().drop(["index"], axis =1)

value_counts = dfFilt['ConsTest'].value_counts()
print(value_counts)

# X represents the features, and y represents the target variable
colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos']
colOpt2 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos']
colOpt3 = ['Var_Al_RelPos', 'Ref_Al_RelPos']
selVars = ['ALLELE.FREQUENCY', "QUAL", "Var_Al_RelPos", "Mean_depth"]

selected_columns = ['Var_Al_RelPos', 'Ref_Al_RelPos']

df_cleaned = dfFilt.dropna(subset=selected_columns).reset_index().drop(["index"], axis = 1)

X = df_cleaned[selVars]
y = df_cleaned['ConsTest']

varNa = X[X.isna().any(axis=1)]
respNa = y[y.isna()]

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Create and train the logistic regression model
model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")

dfFilt.to_csv("snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97_v2.csv", index = False)
