from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import os
from sklearn.metrics import confusion_matrix, roc_auc_score
import numpy as np 

os.chdir("/home/ubuntu/extraVol/ARVAR/iSNVs")

df = pd.read_csv('IntraSnv_results/ampseq_ConsTest_freq_1_0.csv')
#df = pd.read_csv('IntraSnv_ampseq_overlap/ampseq_metaseq_overlap_97_allFreq.csv')
#df = pd.read_csv('IntraSnv_results/metaseq_ConsTest_freq_1_0.csv')

mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)




# make frequency filtering then apply is.in to filter out samples that do not overlapcolOpt1 = ["ALLELE.FREQUENCY", "STRAND.BIAS" , "QUAL", "Var_Al_RelPos", "Ref_Al_RelPos", "meandepth", "coverage",  "meanmapq", "meanbaseq"]
#colOpt1 = ["ALLELE.FREQUENCY", "STRAND.BIAS" , "QUAL", "meandepth", "coverage",  "meanmapq", "meanbaseq"]
colOpt1 = ['ALLELE.FREQUENCY', 'DEPTH', 'Var_Al_RelPos', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq', 'STRAND.BIAS', 'QUAL']
#colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq']


X = dfFilt[colOpt1]
y = dfFilt['ConsTest']

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Combine the features and labels for subsetting
train_data = pd.concat([X_train, pd.Series(y_train, name="Target")], axis=1)

# Separate the dataset into class 0 and class 1
class_0 = train_data[train_data["Target"] == 0]
class_1 = train_data[train_data["Target"] == 1]

# Sample an equal number of samples from both classes
n_samples = min(len(class_0), len(class_1))
class_0_subset = class_0.sample(n=n_samples, random_state=42)
class_1_subset = class_1.sample(n=n_samples, random_state=42)

# Concatenate the subsets to get a balanced training dataset
balanced_train_data = pd.concat([class_0_subset, class_1_subset])

# Split the balanced dataset back into features and labels
X_train_balanced = balanced_train_data.drop("Target", axis=1)
y_train_balanced = balanced_train_data["Target"]

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
model.fit(X_train_balanced, y_train_balanced)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")


print("%.4f" % model.oob_score_)



y_pred = model.predict(X_test)
cm = confusion_matrix(y_test, y_pred)


# Extract the true negative, false positive, false negative, and true positive values
tn, fp, fn, tp = cm.ravel()

# Calculate class error rates
class_0_error = fp / (tn + fp)  # Class 0 (ConsTest = 0)
class_1_error = fn / (fn + tp)  # Class 1 (ConsTest = 1)

# Print the class error rates
print(f"Class 0 Error Rate: {class_0_error:.4f}")
print(f"Class 1 Error Rate: {class_1_error:.4f}")

thresholds = np.linspace(0, 1, 101) 

from sklearn.metrics import precision_score, recall_score, f1_score, average_precision_score

# Calculate model predictions and get predicted probabilities
y_pred_proba = model.predict_proba(X_test)[:, 1]
thresholds = np.linspace(0, 1, 101)  # Range from 0 to 1 with 0.01 increments

results = []

for threshold in thresholds:
    y_pred = (y_pred_proba > threshold).astype(int)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    aupr = average_precision_score(y_test, y_pred_proba)

    results.append({
        "Threshold": threshold,
        "Precision": precision,
        "Recall": recall,
        "F1-Score": f1,
        "AUC-PR": aupr,
    })

# Find the best threshold for F1-Score
best_f1_threshold = max(results, key=lambda x: x["F1-Score"])

# Find the best threshold for AUC-PR
best_aucpr_threshold = max(results, key=lambda x: x["AUC-PR"])

# Find the best threshold for Recall
best_recall_threshold = max(results, key=lambda x: x["Recall"])

print("Best F1-Score Threshold:", best_f1_threshold)
print("Best AUC-PR Threshold:", best_aucpr_threshold)
print("Best Recall Threshold:", best_recall_threshold)

from sklearn.metrics import roc_curve, auc
# Calculate the AUC score
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
roc_auc = auc(fpr, tpr)

# Find the threshold that maximizes the Youden's Index (sum of sensitivity and specificity)
optimal_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[optimal_idx]

y_cust_pred = (y_pred_proba > optimal_threshold).astype(int)


