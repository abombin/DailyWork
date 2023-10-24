### main program
import sys
import subprocess
import pandas as pd
import argparse
import math
import os


parser=argparse.ArgumentParser(description='format and filter Lofreq output')

# add arguments
parser.add_argument(
  "-i",
  "--input_path",
  type=str,
  help="Path to the input file, required"
)

parser.add_argument(
  "-r",
  "--reference",
  type=str,
  help="path to a reference table"
)

parser.add_argument(
  "-o",
  "--output_path",
  type=str,
  nargs="?",
  default="filtered.csv",
  help="path to output file"
)

parser.add_argument(
  "-f",
  "--filter_res",
  type=float,
  nargs="?",
  default=None,
  help="filter resullts to include only relative variant position higher or equal to Float"
)

parser.add_argument(
  "-c",
  "--consensus",
  type=float,
  nargs="?",
  default=0.5,
  help="variant frequencies higher than a value will be considered consensus variants, Float"
)

parser.add_argument(
  "-p",
  "--relative_pos",
  type=float,
  nargs="?",
  default=0.0,
  help="relative position threshold for adjusted frequency, Float"
)

parser.add_argument(
  "-s",
  "--split",
  type=str,
  nargs="?",
  default='F',
  help="write variants that failed the position test separately, string T or F"
)

parser.add_argument(
  "-a",
  "--annot_dir",
  type=str,
  nargs="?",
  default='.',
  help="path to the directory with annotation files, default '.' "
)

# parse arguments
args = parser.parse_args()

df = args.input_path
ref = args.reference
output = args.output_path
filt = args.filter_res
cons = args.consensus
rel = args.relative_pos
splitVar = args.split
annot_dir = args.annot_dir
## main function

def getCorPos(resDf, i, refAl, varAl):
  if len(refAl) > len(varAl):
    alPos = resDf.loc[i, 'POSITION'] + 1
  else:
    alPos = resDf.loc[i, 'POSITION']
  return alPos

# handle deletion
def getRelPosDel(refAl, varAl, i, refSub):
  refRefAl = refAl[1] + '-POS'
  if refRefAl == 'N-POS':
    Ref_Al_RelPos = 'NaN'
  else:
    Ref_Al_RelPos = refSub.loc[0, refRefAl]
  # find which Indel column has deletion amd record relative position
  refVarAl = '-'  + refAl[1:]
  if refSub.loc[0, 'Indel1'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel1-POS']
  elif refSub.loc[0, 'Indel2'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel2-POS']
  elif refSub.loc[0, 'Indel3'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel3-POS']
  else:
    Var_Al_Relpos = 'NaN'
  return Ref_Al_RelPos, Var_Al_Relpos

# process insertions
def getRelPosIns(refAl, varAl, i, refSub):
  refRefAl = refAl + '-POS' 
  Ref_Al_RelPos = refSub.loc[0, refRefAl]
  # find which Indel column has insertion and record relative position
  refVarAl = '+'  + varAl[1:]
  if refSub.loc[0, 'Indel1'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel1-POS']
  elif refSub.loc[0, 'Indel2'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel2-POS']
  elif refSub.loc[0, 'Indel3'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel3-POS']
  else:
    Var_Al_Relpos = 'NaN'
  return Ref_Al_RelPos, Var_Al_Relpos

# process substitutions
def getRelPosSub(refAl, varAl, i, refSub):
  refRefAl = refAl + '-POS'
  Ref_Al_RelPos = refSub.loc[0, refRefAl]
  refVarAl = varAl + '-POS'
  Var_Al_Relpos = refSub.loc[0, refVarAl]
  return Ref_Al_RelPos, Var_Al_Relpos

# get corrected and relative positions for all SNVs
def getCorAndRelPos(refDf, resDf):
  corPos = []
  Ref_Al_RelPos = []
  Var_Al_RelPos = []
  for i in range(len(resDf.index)):
    refAl = resDf.loc[i, 'REF-NT']
    varAl = resDf.loc[i, 'VAR-NT']
    curCorPos = getCorPos(resDf = resDf, i = i, refAl = refAl, varAl = varAl)
    refSub = refDf.loc[refDf['POS'] == curCorPos].reset_index()
    if len(refAl) > len(varAl):
      try:
        curRefPos, curVarPos = getRelPosDel(refAl = refAl, varAl= varAl, i = i, refSub = refSub)
      except:
        print(f'{i} result table index deletion failed')
        sys.exit(1)
    elif len(refAl) < len(varAl):
      try:
        curRefPos, curVarPos = getRelPosIns(refAl = refAl, varAl= varAl, i = i, refSub = refSub)
      except:
        print(f'{i} result table index insertion failed')
        sys.exit(1)
    elif len(refAl) == len(varAl):
      try:
        curRefPos, curVarPos = getRelPosSub(refAl = refAl, varAl= varAl, i = i, refSub = refSub)
      except:
        print(f'{i} result table index substitution failed')
        sys.exit(1)
    corPos.append(curCorPos)
    Ref_Al_RelPos.append(curRefPos)
    Var_Al_RelPos.append(curVarPos)
  resDf['Position_corrected'] = corPos
  resDf['Ref_Al_RelPos'] = Ref_Al_RelPos
  resDf['Var_Al_RelPos'] = Var_Al_RelPos
  return resDf

# write results
def writeResults(resDf, splitVar, output):
  if splitVar == 'T':
    posRes = finalDf.loc[finalDf['Position_test'] == 'Pass']
    negRes = finalDf.loc[finalDf['Position_test'] == 'Fail' ]
    negRes.to_csv(path_or_buf= 'failed_position_test.csv', index=False)
    posRes.to_csv(path_or_buf= output, index=False)
  else:
    resDf.to_csv(path_or_buf= output, index=False)

# runn all functions
def runAll():
  current_directory = os.getcwd()
  try:
    refDf=pd.read_table(ref, sep = '\t')
    resDf=pd.read_table(df, sep = '\t')
  except:
    print(f"Failed to read files {current_directory}")
    sys.exit(1)
  try:
    resDf = getCorAndRelPos(refDf = refDf, resDf = resDf)
  except:
    print(f"Failed to find corrected and relative positions {current_directory}")
    sys.exit(1)
  return resDf

# execute runAll function
if __name__ == "__main__":
  resDf = runAll()
  writeResults(resDf = resDf, splitVar = splitVar, output = output)
