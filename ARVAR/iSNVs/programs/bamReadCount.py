#lofreq 2.1.5
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
  "-m",
  "--mapqual",
  type=str,
  nargs="?",
  default='20',
  help="Mapping quality for bam-readcount, default 20"
)

parser.add_argument(
  "-q",
  "--basequal",
  type=str,
   nargs="?",
   default='20',
  help="Base quality for bam-readcount, default 20"
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
mapqual = args.mapqual
basequal = args.basequal
lofreq_filter = args.allele_frequency
prefix = args.prefix
bin_path = args.bin_path

startTime = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
roseColNames = ['POS', 'A-POS', 'T-POS', 'C-POS', 'G-POS', 'Indel1', 'Indel1-POS', 'Indel2', 'Indel2-POS', 'Indel3', 'Indel3-POS']
curHeaders = ['position', 'A_avg_pos_as_fraction', 'T_avg_pos_as_fraction', 'C_avg_pos_as_fraction', 'G_avg_pos_as_fraction', \
              'Indel1_base', 'Indel1_avg_pos_as_fraction', 'Indel2_base', 'Indel2_avg_pos_as_fraction', 'Indel3_base', 'Indel3_avg_pos_as_fraction']
              
# get reference name and length from alignment file bowtie2
def getRefLen(f):
  with open(f) as input_file:
    head = list(islice(input_file, 2))
  targLine = head[1].split('\t')
  ref = targLine[1].replace("SN:", "")
  alLen = targLine[2].replace("LN:", "").replace("\n", "")
  return ref, alLen



# bam readcount
def bamReadCount(mapqual, basequal, refFasta, bam, refName, refLen, bin_path):
  cmd_str = '%s/bam-readcount -w1 -q %s -b %s -f %s %s %s:1-%s > readcount.tsv' % (bin_path, mapqual, basequal, refFasta, bam, refName, refLen)
  subprocess.run(cmd_str, shell=True)

# parse readcount file 
def getReadCountRes():
  with open("readcount.tsv", "r") as f:
    lines=f.readlines()
    result=[]
    for x in lines:
      result.append(x.split('\t'))
  return(result)

# get number of columns/indels
def getColNumbers(readCountRes):
  colNumb =[]
  for i in readCountRes:
    colNumb.append(len(i))
  uniqNumb = list(np.unique(colNumb))
  return uniqNumb

# extract info for nucleotides
def parseNucleotide(i, col, readCountRes):
  nuclRes = readCountRes[i][col].split(':')
  selRes = list(nuclRes[1].split(' ')) + nuclRes[5:9]
  return selRes

# extract infor for indels
def parseIndels(i, col, readCountRes):
  emptList = ["NaN"] * 6
  if len(readCountRes[i]) > col:
    nuclRes = readCountRes[i][col].split(':')
    selRes = nuclRes[0:2] + nuclRes[5:9]
  else:
    selRes = emptList
  return selRes

# parse whole bam readcount
def parseReadCount(readCountRes):
  allMainRes = []
  allAres = []
  allCres = []
  allGres = []
  allTres = []
  allIndel1 = []
  allIndel2 = []
  allIndel3 = []
  for i in range(len(readCountRes)):
    mainRes = readCountRes[i][0:4]
    aRes = parseNucleotide(i = i, col = 5, readCountRes = readCountRes)
    cRes = parseNucleotide(i = i, col = 6, readCountRes = readCountRes)
    gRes = parseNucleotide(i = i, col = 7, readCountRes = readCountRes)
    tRes = parseNucleotide(i = i, col = 8, readCountRes = readCountRes)
    indel1 = parseIndels(i = i , col = 10, readCountRes = readCountRes)
    indel2 = parseIndels(i = i , col = 11, readCountRes = readCountRes)
    indel3 = parseIndels(i = i , col = 12, readCountRes = readCountRes)
    #combine lists
    allMainRes.append(mainRes)
    allAres.append(aRes)
    allCres.append(cRes)
    allGres.append(gRes)
    allTres.append(tRes)
    allIndel1.append(indel1)
    allIndel2.append(indel2)
    allIndel3.append(indel3)
    
  # make dataframes  
  mainData = pd.DataFrame(allMainRes, columns = ['chr',	'position',	'reference_base',	'depth'])
  aDf = pd.DataFrame(allAres, columns = ['A_count', 'A_num_plus_strand', 'A_num_minus_strand', 'A_avg_pos_as_fraction', 'A_avg_num_mismatches_as_fraction'])
  cDf = pd.DataFrame(allCres, columns = ['C_count', 'C_num_plus_strand', 'C_num_minus_strand', 'C_avg_pos_as_fraction', 'C_avg_num_mismatches_as_fraction'])
  gDf = pd.DataFrame(allGres, columns = ['G_count', 'G_num_plus_strand', 'G_num_minus_strand', 'G_avg_pos_as_fraction', 'G_avg_num_mismatches_as_fraction'])
  tDf = pd.DataFrame(allTres, columns = ['T_count', 'T_num_plus_strand', 'T_num_minus_strand', 'T_avg_pos_as_fraction', 'T_avg_num_mismatches_as_fraction'])
  indel1Df = pd.DataFrame(allIndel1, columns = ['Indel1_base', 'Indel1_count', 'Indel1_num_plus_strand', 'Indel1_num_minus_strand', 'Indel1_avg_pos_as_fraction', 'Indel1_avg_num_mismatches_as_fraction'])
  indel2Df = pd.DataFrame(allIndel2, columns = ['Indel2_base', 'Indel2_count', 'Indel2_num_plus_strand', 'Indel2_num_minus_strand', 'Indel2_avg_pos_as_fraction', 'Indel2_avg_num_mismatches_as_fraction'])
  indel3Df = pd.DataFrame(allIndel3, columns = ['Indel3_base', 'Indel3_count', 'Indel3_num_plus_strand', 'Indel3_num_minus_strand', 'Indel3_avg_pos_as_fraction', 'Indel3_avg_num_mismatches_as_fraction'])
  # combine tables  
  combData = pd.concat([mainData, aDf, cDf, gDf, tDf, indel1Df, indel2Df, indel3Df], axis="columns")
  return combData

# write tables
def writeDf(df, outName):
  if len(df.index) > 0:
    df.to_csv(path_or_buf = outName, index = False, sep = "\t")
    
# run the pipeline  
def runAll():
  #

  refName, refLen = "MN908947.3", "29903"

  #

  try:
    bamReadCount(mapqual = mapqual, basequal = basequal, refFasta = refFasta, bam = 'output.bam', refName = refName, refLen = refLen, bin_path = bin_path)
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    readCountRes = getReadCountRes()
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    mainData = parseReadCount(readCountRes = readCountRes)
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    writeDf(df = mainData, outName = "readCountParsed.tsv")
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    selectData = mainData[curHeaders]
    selectData.columns = roseColNames
    writeDf(df = selectData, outName = prefix + "_pos-filter.tsv")
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  print("DONE!")

# run tool
if __name__ == "__main__":
  runAll()
    
