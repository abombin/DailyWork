import subprocess
import os
import pandas as pd

os.chdir("/home/ubuntu/extraVol/ARVAR/iSNVs")

os.makedirs("dx_metaseq_paths", exist_ok=True)

# get samples list
def getSamplesNames(df, curType):
  allSampleNames = []
  curDf = pd.read_csv(df)
  for i in range(len(curDf.index)):
    if curType == "metaseq":
      curSeq = curDf.loc[i, "MetaseqNames"]
    if curSeq  != "NA" and not pd.isna(curSeq):
      curSample = curDf.loc[i, "Sample_id"]
      allSampleNames.append(curSample)
  return allSampleNames

# download paths for selected samples
def writeSamplesDx(curSamples, outDir):
  for curSample in curSamples:
    curName = curSample.replace("-", "*").replace("_", "*")
    searchName = "*" + curName + "*"
    cmd_str = f'dx find data --name "{searchName}"'
    result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
    output = result.stdout
    outFile = f'{outDir}/{curSample}'
    with open(outFile, 'w') as file:
      file.write(output)


curSamples = getSamplesNames(df = "SeqSamples_Eval_Depth10.csv", curType = "metaseq")
writeSamplesDx(curSamples = curSamples, outDir = "dx_metaseq_paths")   
