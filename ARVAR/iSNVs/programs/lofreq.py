import argparse
import os
import subprocess
from itertools import islice
import pandas as pd
import numpy as np
import glob
from datetime import datetime
import sys
import re
import shutil

parser=argparse.ArgumentParser(description='get and parse bam-readcount and lofreq outputs')

# add arguments

parser.add_argument(
  "-r",
  "--reference_path",
  type=str,
  help="Path to the reference file, required",
  metavar=""
)

parser.add_argument(
  "-p",
  "--prefix",
  type=str,
  nargs="?",
  default='sample',
  help="Prefix for the output files, default sample"
)


parser.add_argument(
  "-a",
  "--allele_frequency",
  type=float,
   nargs="?",
   default=0.02,
  help="Minimal allele frequency to keep the variant in the final output, default 0.02"
)

parser.add_argument(
  "-t",
  "--threads",
  type=str,
   nargs="?",
   default='8',
  help="Threads to use, default 8"
)

parser.add_argument(
  "-b",
  "--bin_path",
  type=str,
   nargs="?",
   default='programs/bin',
  help="Path to bam-readcount and lofreq binaries"
)

# parse arguments
args = parser.parse_args()

refFasta = args.reference_path
threads = args.threads
lofreq_filter = args.allele_frequency
prefix = args.prefix
bin_path = args.bin_path

# write tables
def writeDf(df, outName):
  if len(df.index) > 0:
    df.to_csv(path_or_buf = outName, index = False, sep = "\t")
  
# prep for lowfreq
def indelqual(bam, refFasta, bin_path):
  cmd_str = f"{bin_path}/lofreq indelqual --dindel -f {refFasta} -o sample_dindel.tmp.bam {bam}"
  subprocess.run(cmd_str, shell = True)
  sam_str = f"{bin_path}/samtools index sample_dindel.tmp.bam"
  subprocess.run(sam_str, shell = True)

# call lowfreq
def lofreq(bam, refFasta, t, bin_path):
  cmd_str = f"{bin_path}/lofreq call --call-indels -f {refFasta} -o sample_lf.vcf {bam} 2>sample_lofreq-log.txt"
  subprocess.run(cmd_str, shell = True)

# read lowfreq results
def getlowfreqRes():
  with open("sample_lf.vcf", "r", errors='ignore') as f:
    lines=f.readlines()
    result=[]
    for i in range(len(lines)):
      curLine = lines[i].replace("\n", "")
      if curLine[0] != "#":
        curLineList = re.split('\t|;|=|,', curLine)
        result.append(curLineList)
  return(result)

# subset list by index
def subsetList(indList, curList):
  subList = []
  for i in indList:
    subList.append(curList[i])
  return subList

# filter by the number of reference reads
def filterRefNumb(lofreqDf, sumReads):
  lofreqDf["FWD-REF"] = lofreqDf["FWD-REF"].astype(int)
  lofreqDf["REV-REF"] = lofreqDf["REV-REF"].astype(int)
  dfFilter = lofreqDf.loc[(lofreqDf["FWD-REF"] + lofreqDf["REV-REF"]) >= sumReads]
  return(dfFilter)

# parse lowfreq
def parseLowfreq(lofreqRes):
  if len(lofreqRes) == 0:
    print("No Variants Were identified")
  colNames = ['REF-GENOME','POSITION','REF-NT','VAR-NT','QUAL','FILTER','DEPTH','ALLELE-FREQUENCY','STRAND-BIAS','FWD-REF','REV-REF','FWD-VAR','REV-VAR']
  indList = [0,1,3,4,5,6,8,10,12,14,15,16,17]
  results = []
  for i in range(len(lofreqRes)):
    try:
      curLine = lofreqRes[i]
      #if curLine[6] == "PASS" and float(curLine[10]) >= lofreq_filter:  # for now remove pass parameter since some snps are missing
      if len(curLine) > 9 and float(curLine[10]) >= lofreq_filter:
        subsRes = subsetList(indList = indList, curList = curLine)
        results.append(subsRes)
    except Exception as e:
      print(e)
  lofreqDf = pd.DataFrame(results, columns = colNames)
  lofreqDf = filterRefNumb(lofreqDf, 0) # replaced 2 with 0 to find missing snps
  return lofreqDf

# clean directory
def cleanDir():
  os.makedirs("process_files/", exist_ok = True)
  delList = glob.glob("*dindel*")
  for i in delList:
    os.remove(i)
  mv_list = glob.glob("*bam*")
  mv_list.extend(["readcount.tsv", "sample_lofreq-log.txt", "sample_lf.vcf", "readCountParsed.tsv", "editFile.sam"])
  for i in range(len(mv_list)):
    curFile = str(mv_list[i])
    if os.path.isfile(curFile):
      shutil.move(curFile, "process_files/")
  
def writeSummary(startTime, mapqual):
  curDate = datetime.today().strftime('%Y-%m-%d')
  endTime = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
  with open("runSummary.txt", "a") as text_file:
    text_file.write(f"Start Date {startTime}" + "\n")
    text_file.write(f"End Date {endTime}" + "\n")
    text_file.write(f"Mapping quality {mapqual}" + "\n")

# if file exists delete it
def delFile(fPath):
  if os.path.isfile(fPath):
    os.remove(fPath)



def runAll():
  #
  delFile(fPath = "sample_dindel.tmp.bam")
  delFile(fPath = "sample_dindel.tmp.bam.bai")
  #
  try:
    indelqual(bam = 'output.bam', refFasta = refFasta, bin_path = bin_path)
  except Exception as e:
    print(e)
    print("Lofreq indels quality annotations failed")
    sys.exit(1)
  #
  try:
    lofreq(bam = 'sample_dindel.tmp.bam', refFasta = refFasta, t = threads, bin_path = bin_path)
  except Exception as e:
    print(e)
    print("Lofreq run failed")
    sys.exit(1)
  #
  if os.path.exists("sample_lf.vcf"):
    pass
  else:
    lofreq(bam = 'output.bam', refFasta = refFasta, t = threads, bin_path = bin_path)
  #
  try:
    lofreqRes = getlowfreqRes()
  except Exception as e:
    print(e)
    print("Lofreq file parsing failed")
    sys.exit(1)
  #
  try:
    lowfreqDf = parseLowfreq(lofreqRes = lofreqRes)
    writeDf(df = lowfreqDf, outName = prefix + "_lofreq-output.tsv")
  except Exception as e:
    print(e)
    print("Lofreq file annotation failed")
    sys.exit(1)
  # wrap up
  #cleanDir()
  #writeSummary(startTime = startTime, mapqual = mapqual)
  print("DONE!")

# run tool
if __name__ == "__main__":
  runAll()
