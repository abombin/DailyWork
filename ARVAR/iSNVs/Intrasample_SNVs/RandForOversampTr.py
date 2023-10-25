from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, roc_auc_score
from imblearn.over_sampling import SMOTE
import pandas as pd

# Load your dataset
df = pd.read_csv('IntraSnv_results/metaseq_ConsTest_freq_1_0.csv')

mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)

# make frequency filtering then apply is.in to filter out samples that do not overlapcolOpt1 = ["ALLELE.FREQUENCY", "STRAND.BIAS" , "QUAL", "Var_Al_RelPos", "Ref_Al_RelPos", "meandepth", "coverage",  "meanmapq", "meanbaseq"]
colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq']
colOpt5 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos',  'meandepth', 'meanbaseq', 'Sample', 'Sample_AlignPos_Ref_Var']

X = dfFilt[colOpt5]
y = dfFilt['ConsTest']

varNa = X[X.isna().any(axis=1)]
respNa = y[y.isna()]

X_train_1, X_test_1, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

X_train = X_train_1[colOpt1]
X_test = X_test_1[colOpt1]

# Apply SMOTE to the training dataset
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

# Create and fit the Random Forest model on the resampled training dataset
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
model.fit(X_train_resampled, y_train_resampled)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]
y_pred = model.predict(X_test)

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")

cm = confusion_matrix(y_test, y_pred)


# Extract the true negative, false positive, false negative, and true positive values
tn, fp, fn, tp = cm.ravel()

# Calculate class error rates
class_0_error = fp / (tn + fp)  # Class 0 (ConsTest = 0)
class_1_error = fn / (fn + tp)  # Class 1 (ConsTest = 1)

# Print the class error rates
print(f"Class 0 Error Rate: {class_0_error:.4f}")
print(f"Class 1 Error Rate: {class_1_error:.4f}")

print(f"AUC Score: {auc_score}")