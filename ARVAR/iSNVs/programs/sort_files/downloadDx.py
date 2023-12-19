import os
import pandas as pd
import subprocess
from multiprocessing import Pool

os.makedirs("metaseq_pass_fastqs", exist_ok = True)

passDf = pd.read_csv("metaseq_libs_pass.csv")
passList = passDf["combFastq"].to_list()

resolvDf  = pd.read_csv("meta_lib_reolv_spike_size.csv")
resolvList = resolvDf["V6"].to_list()

def splitFiles(curList):
  combList = []
  for i in curList:
    curElem = i.split(";")
    combList.extend(curElem)
  return(combList)

passList = splitFiles(curList = passList)

combList = passList + resolvList

def downloadFiles(inFile):
  cmd_str = f"dx download {inFile} -o metaseq_pass_fastqs/"
  subprocess.run(cmd_str, shell = True)

with Pool(7) as pool:
  pool.map(downloadFiles, combList)
