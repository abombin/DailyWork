import os
import subprocess
import glob
import shutil
from multiprocessing import Pool
from functools import partial

t = 1
tSmall = 1
suffix = "_001"

# currently use 4 cores for multiprocessor parallel tasks

index_ref = "false"


os.makedirs("preprocess", exist_ok = True)
os.makedirs("input", exist_ok = True)

# index reference
def indRef(index_ref):
  if index_ref == "true":
    os.chdir("references")
    cmd_str = f"hisat2-build -p 16 SARSCov2_Ref.fasta SARSCov2_Ref"
    subprocess.run(cmd_str, shell = True)
    os.chdir("../")

indRef(index_ref = index_ref)
  
# 
def getSamplesName(suffix):
  filesList = glob.glob("./fastqs/*_R1*")
  samplesName = []
  for i in filesList:
    newName = i.replace("./fastqs/", "").replace(f'_R1{suffix}.fastq.gz', "")
    samplesName.append(newName)
  return samplesName

samplesNames = getSamplesName(suffix = suffix)

def splitFiles(sampleName):
  targDir = f'preprocess/{sampleName}/'
  os.makedirs(targDir, exist_ok = True)
  curSamples = glob.glob(f"fastqs/{sampleName}*")
  for curSample in curSamples:
      shutil.copy(curSample, targDir)

with Pool(tSmall) as pool:
  pool.map(splitFiles, samplesNames)

print("Done copying files!")

# filter fastqs with fastp seems like 2 threads are not enough
def runFastp(suffix, sampleName):
  targDir = "preprocess/" + sampleName + "/"
  os.chdir(targDir)
  cmd_str = f'../../programs/bin/fastp -i {sampleName}_R1{suffix}.fastq.gz -I {sampleName}_R2{suffix}.fastq.gz \
  -o filtered_1.fastq -O filtered_2.fastq -D -A -l 25 -w 4 -h fastp.html'
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

with Pool(t) as pool:
  suffix = suffix
  fastp = partial(runFastp, suffix)
  pool.map(fastp, samplesNames)

print("Done fastp files!")  

def runHisat2(sampleName):
  targDir = "preprocess/" + sampleName + "/"
  os.chdir(targDir)
  cmd_str = f"hisat2 -p 4 -x ../../references/SARSCov2_Ref -1 filtered_1.fastq -2 filtered_2.fastq -S {sampleName}.sam"
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

with Pool(t) as pool:
  pool.map(runHisat2, samplesNames)

print("Done Hisat2 files!")  

def copyFiles(sampleName):
  targDir = "preprocess/" + sampleName + "/"
  os.chdir(targDir)
  shutil.move(f'{sampleName}.sam', "../../input/")
  os.chdir("../../")

with Pool(tSmall) as pool:
  pool.map(copyFiles, samplesNames)

print("Preprocessing complete!")
