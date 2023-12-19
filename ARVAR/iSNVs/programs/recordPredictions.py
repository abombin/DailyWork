import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

# Define the file path
DF = "snvs_comb_res/ampseq_metaseq_overlap_comb_derep_decont_covFilt_97.csv"

# Set the random seed
curSeed = 42

# Read the CSV file into a DataFrame
df = pd.read_csv(DF)

# Filter out rows with NA values in Var_Al_RelPos
dfFilt = df.dropna(subset=['Var_Al_RelPos'])
dfFilt['Var_Al_RelPos'] = pd.to_numeric(dfFilt['Var_Al_RelPos'], errors='coerce')
dfFilt['Ref_Al_RelPos'] = pd.to_numeric(dfFilt['Ref_Al_RelPos'], errors='coerce')

# Set the random seed
np.random.seed(curSeed)

# Split the data into training and test sets
train_data, test_data = train_test_split(dfFilt, test_size=0.3, random_state=curSeed)

# Fit the logistic regression model
features = ['ALLELE.FREQUENCY', 'DEPTH', 'Var_Al_RelPos', 'Mean_depth']
target = 'ConsTest'

X_train = train_data[features]
y_train = train_data[target]

X_test = test_data[features]

model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

# Predict probabilities on the test set
probs = model.predict_proba(X_test)[:, 1]

# Define the threshold
threshold = 0.5

# Create a new column for predicted values based on the threshold
test_data['PredictedConsTest'] = np.where(probs >= threshold, 1, 0)

# Calculate the number of correct predictions
equal_count = sum(test_data['ConsTest'] == test_data['PredictedConsTest'])

# Calculate the number of incorrect predictions
not_equal_count = sum(test_data['ConsTest'] != test_data['PredictedConsTest'])

# Calculate the accuracy
sumCorrect = equal_count / (equal_count + not_equal_count) * 100
print("Accuracy:", sumCorrect)

# Calculate the AUC value
roc_auc = roc_auc_score(test_data['ConsTest'], probs)
print("AUC:", roc_auc)
