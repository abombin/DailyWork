import pandas as pd
import glob
import shutil
import os

suffix = "_001"

# get names of the sampels from the input folder
def getSamplesName(inDir,suffix):
  filesList = glob.glob(f"./{inDir}/*_R1*")
  samplesName = []
  for i in filesList:
    newName = i.replace(f"./{inDir}/", "").replace(f'_R1{suffix}.fastq.gz', "")
    samplesName.append(newName)
  return samplesName

# make input csv file
def inputDf(inDir, samplesList):
  r1 = []
  r2 = []
  for curSample in samplesList:
    curR1 = glob.glob(f"{inDir}/{curSample}*R1*")[0]
    curR2 = glob.glob(f"{inDir}/{curSample}*R2*")[0]
    r1.append(curR1)
    r2.append(curR2)
  combDf = pd.DataFrame({"sample":samplesList, "fastq_1":r1, "fastq_2":r2})
  return combDf

def writeInputDat(df, outName):
  if len(pd.unique(df['sample'])) == len(df.index):
    df.to_csv(path_or_buf = outName, index = False)
    print('input recorded')

# regular samples
samplesNames = getSamplesName(inDir = "fastqs", suffix = suffix)
inputDat = inputDf(inDir = "fastqs", samplesList = samplesNames)
writeInputDat(df=inputDat, outName="samples_input.csv")

# water samples
waterSamples = getSamplesName(inDir = "water_fastqs", suffix = suffix)
waterDat = inputDf(inDir = "water_fastqs", samplesList = waterSamples)
writeInputDat(df=waterDat, outName="water_input.csv")

# pd.set_option('display.max_colwidth', None)