import subprocess
import os
import glob
import shutil
from multiprocessing import Pool
import time

t = 2

## get iSNV files
filesList = glob.glob("input/*.sam")

def prepFiles():
  samplesList = []
  filesList = glob.glob("input/*.sam")
  for file in filesList:
    sampleName = file.replace("input/", "").replace(".sam", "")
    targDir = "process/" + sampleName + "/"
    samplesList.append(sampleName)
    os.makedirs(targDir, exist_ok = True)
    shutil.copy(file, targDir)
  return samplesList

samplesList = prepFiles()


def runGetVarFiles(sampleName):
  targDir = "process/" + sampleName + "/"
  os.chdir(targDir)
  sample = sampleName + ".sam"
  cmd_str = f'python ../../programs/src/getVarFilesCLT.py -i {sample} \
-r ../../programs/data/SARSCov2_Ref.fasta -b ../../programs/bin -t 4'
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

# run lofreq pipe
start_time = time.time()
with Pool(t) as pool:
  pool.map(runGetVarFiles, samplesList)

print("getVarFiles completed in --- %s minutes ---" % ((time.time() - start_time) /60) )




  
