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
  "-i", "--input_path",
  type=str,
  metavar='',
  help="Path to the input file, required"
)

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
  default='0',
  help="Mapping quality for bam-readcount, default 0"
)

parser.add_argument(
  "-q",
  "--basequal",
  type=str,
   nargs="?",
   default='0',
  help="Base quality for bam-readcount, default 0"
)

parser.add_argument(
  "-a",
  "--allele_frequency",
  type=float,
   nargs="?",
   default=0.0,
  help="Minimal allele frequency to keep the variant in the final output, default 0"
)

parser.add_argument(
  "-t",
  "--threads",
  type=str,
   nargs="?",
   default='1',
  help="Threads to use, default 1"
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

inFile = args.input_path
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
              
# main funtions 

# get reference name and length from alignment file bowtie2
def getRefLen(f):
  with open(f) as file:
    for line in file:
      curLine = line.split('\t')
      if curLine[0] == "@SQ":
        ref = curLine[1].replace("SN:", "")
        alLen = curLine[2].replace("LN:", "").replace("\n", "")
        break
  return ref, alLen

# replcase N with D in cigar string
def editSam(inFile):
  with open(inFile) as file, open('editFile.sam','w') as secondfile:
      for line in file:
          curLine = line.split('\t')
          if curLine[0][0] == "@":
            newline = "\t".join(curLine)
          elif len(curLine) < 5:
            newline = "\t".join(curLine)
          elif len(curLine) > 5:
            curLine[5] = curLine[5].replace("N", "D")
            newline = "\t".join(curLine)
          secondfile.write(newline)


#Convert the sam file into a bam file containing only mapped reads
def sam2Bam(sam, threads, bin_path):
  cmd_str = f"{bin_path}/samtools view -@ {threads} -bu -F 4 {sam} -o - | {bin_path}/samtools sort -@ {threads} - -o output.bam"
  subprocess.run(["bash", "-c", cmd_str],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          text=True)
                          
# index bam file
def indexBam(bam, threads, bin_path):
  cmd_str = f'{bin_path}/samtools index -@ {threads} {bam}'
  subprocess.run(cmd_str, shell=True, stderr=subprocess.PIPE)

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

# run the pipeline  
def runAll():
  #
  try:
    refName, refLen = getRefLen(f = inFile)
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    editSam(inFile = inFile)
  except Exception as e:
    print(e)
    sys.exit(1)
    
  try:
    sam2Bam(sam = 'editFile.sam', threads = threads, bin_path = bin_path)
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    indexBam(bam = 'output.bam', threads = threads, bin_path = bin_path)
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    bamReadCount(mapqual = mapqual, basequal = basequal, refFasta = refFasta, bam = 'output.bam', refName = refName, refLen = refLen, bin_path = bin_path)
  except Exception as e:
    print(e)
    sys.exit(1)
  #
  try:
    os.remove("output.sam")
    os.remove('editFile.sam')
  except:
    pass
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
  writeSummary(startTime = startTime, mapqual = mapqual)
  print("DONE!")

# run tool
if __name__ == "__main__":
  runAll()
    
