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
  default=0.2,
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

# adjust variants frequency based on threshold
def adjFreq(df, relPos):
  newFreq = []
  allTests = []
  df['Ref_Al_RelPos'] = df['Ref_Al_RelPos'].astype(float)
  df['Var_Al_RelPos'] = df['Var_Al_RelPos'].astype(float)
  for i in range(len(df.index)):
    if df.loc[i, 'Level'] == 'consensus' and df.loc[i, 'Ref_Al_RelPos'] >= relPos:
      freq = df.loc[i, 'ALLELE-FREQUENCY']
      posTest = 'Pass'
    elif df.loc[i, 'Level'] == 'consensus' and df.loc[i, 'Ref_Al_RelPos'] < relPos:
      freq = 1
      posTest = 'Fail'
    elif df.loc[i, 'Level'] == 'minority' and df.loc[i, 'Var_Al_RelPos'] >= relPos:
      freq = df.loc[i, 'ALLELE-FREQUENCY']
      posTest = 'Pass'
    elif df.loc[i, 'Level'] == 'minority' and df.loc[i, 'Var_Al_RelPos'] < relPos:
      freq = 0
      posTest = 'Fail'
    newFreq.append(freq)
    allTests.append(posTest)
  df['Position_test'] = allTests
  df['Freq_adj'] = newFreq
  return(df)

# find if the variant is a consensus level
def isConsensus(df, freq):
  consList = []
  for i in range(len(df.index)):
    if df.loc[i, 'ALLELE-FREQUENCY'] > freq:
      consList.append('consensus')
    else:
      consList.append('minority')
  df['Level'] = consList
  return(df)

# correct reference and variant alleles, and add nucleotide change
def typeAllele(df):
  newType = []
  corRefAl = []
  corVarAl = []
  allNuclCh =[]
  for i in range(len(df.index)):
    refAl = df.loc[i, 'REF-NT']
    varAl = df.loc[i, 'VAR-NT']
    if len(refAl) > len(varAl):
      typeAl = 'Deletion'
      newRefAl = refAl[1] # confirm that reference is always a second character
      newVarAl = '-' + refAl[1:]
      nuclCh = str(df.loc[i, 'Position_corrected']) + newVarAl
    elif len(refAl) < len(varAl):
      typeAl = 'Insertion'
      newRefAl = refAl 
      newVarAl = '+' + varAl[1:]
      nuclCh = str(df.loc[i, 'Position_corrected']) + newVarAl
    elif len(refAl) == len(varAl):
      typeAl = 'Substitution'
      newRefAl = refAl 
      newVarAl = varAl
      nuclCh = refAl + str(df.loc[i, 'Position_corrected']) + varAl
    newType.append(typeAl)
    corRefAl.append(newRefAl)
    corVarAl.append(newVarAl)
    allNuclCh.append(nuclCh)
  df['Type'] = newType
  df['REF-allele_corrected'] = corRefAl
  df['VAR-allele_corrected'] = corVarAl
  df['Nucleotide_Change'] = allNuclCh
  return df

# annotate substitutions
def annotSubst(pos, codNumb, region, corVarAl, annot, codons):
  # subset codon table
  selAnot = annot.loc[(annot['All_Products'] == region) & (annot['Codon#'] == str(codNumb))].reset_index().drop(['index'], axis=1)
  if len(selAnot.index) < 3:
    mType = 'Synonymous'
    chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
    oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
    newAA = oldAA
  else:
    # find index of the nucleotide to change 
    chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
    # change nucleotide
    selAnot.loc[chInd, 'Ref'] = corVarAl
    # write new codon
    codList = selAnot['Ref'].tolist()
    newCod = ''.join(codList)
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
    oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
    if newAA == oldAA:
      mType = 'Synonymous'
    else:
      mType = 'Nonsynonymous'
  aaChange = oldAA + codNumb + newAA
  return [mType, aaChange]

# annotate insertions
def annotInsert(pos, codNumb, region, corVarAl, annot, codons):
  nucleotides = corVarAl[1:]
  # subset codon table
  selAnot = annot.loc[(annot['All_Products'] == region) & (annot['Codon#'] == str(codNumb))].reset_index().drop(['index'], axis=1)
  # find index of the nucleotide to change 
  chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
  # get reference amino acid
  oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
  # find new codon and aminoacid
  if chInd == 2:
    newAA = oldAA
  elif chInd == 1:
    codList = selAnot['Ref'].tolist()
    codList.insert(2, nucleotides)
    newSeq = ''.join(codList)
    newCod = newSeq[0:3]
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  elif chInd == 0:
    codList = selAnot['Ref'].tolist()
    codList.insert(1, nucleotides)
    newSeq = ''.join(codList)
    newCod = newSeq[0:3]
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  # inframe or not  
  if len(nucleotides) % 3 == 0:
    mType = 'Inframe'
    addNumb = '+' + str(int(len(nucleotides)/3)) + 'AA'
  else:
    mType = 'Frameshift'
    addNumb = '+' + str(int(len(nucleotides))) + 'nt'
  # if to write a new aa  
  if chInd == 2:
    aaChange = ''.join([oldAA, codNumb, addNumb])
  else:
    aaChange = ''.join([oldAA, codNumb, newAA, addNumb])
  return [mType, aaChange]

# annotate deletions
def annotDel(pos, codNumb, region, corVarAl, annot):
  nucleotides = corVarAl[1:]
  # subset codon table
  selAnot = annot.loc[(annot['All_Products'] == region) & (annot['Codon#'] == str(codNumb))].reset_index().drop(['index'], axis=1)
  # find index of the nucleotide to change 
  chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
  # deletions never need a new amino acid or a check at what codon position nucleotide is deleted
  # needs old amino acid, codon number, -numb of aa or nt
  oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
  if len(nucleotides) % 3 == 0:
    mType = 'Inframe'
    addNumb = '-' + str(int(len(nucleotides)/3)) + 'AA'
  else:
    mType = 'Frameshift'
    addNumb = '-' + str(int(len(nucleotides))) + 'nt'
  aaChange = ''.join([oldAA, codNumb, addNumb])
  return [mType, aaChange]

# use functions to annotate the table
def annotate(df, annot, codons):
  allRegions = []
  allMutTypes = []
  allAaCh = []
  for i in range(len(df.index)):
    pos = df.loc[i, 'Position_corrected']
    region = annot.loc[annot['NT'] == pos, 'All_Products'].iloc[0]
    codNumb = annot.loc[annot['NT'] == pos, 'Codon#'].iloc[0]
    corVarAl = df.loc[i, 'VAR-allele_corrected']
    newRegion = annot.loc[annot['NT'] == pos, 'All_Products_Edit'].iloc[0]
    if region == "5'UTR" or region == "3'UTR" or region == 'E/NCR' or region == 'NaN' or pd.isna(region) or region == 'NCR':
      mutType = 'Synonymous'
      aaCh = "NCR"
      if newRegion == "NaN" or pd.isna(newRegion):
        newRegion = "NCR"
    else:
      # handle substitutions
      if df.loc[i, 'Type'] == 'Substitution':
        resList = annotSubst(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl, annot = annot, codons = codons)
      # handle insertions  
      elif df.loc[i, 'Type'] == 'Insertion':
        resList = annotInsert(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl, annot = annot, codons = codons)
      # handle deletions
      elif df.loc[i, 'Type'] == 'Deletion':
        resList = annotDel(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl, annot = annot)
      mutType = resList[0]
      aaCh = resList[1]
    # combine lists
    allRegions.append(newRegion)
    allMutTypes.append(mutType)
    allAaCh.append(aaCh)
    # make new columns
  df['Region'] = allRegions
  df['Mutation_type'] = allMutTypes
  df['AA_change'] = allAaCh
  return df

# add  Pi*Ln(Pi) for shannon index 
def ShPi(df):
  pi = []
  for i in range(len(df.index)):
    if df.loc[i, 'Freq_adj'] == 0:
      newPi = 0
    else:
      newPi = df.loc[i, 'Freq_adj'] * math.log(df.loc[i, 'Freq_adj'])
    pi.append(newPi)
  df['Pi*Ln(Pi)'] = pi
  return df

# filter results based on relative position
def filterResRelPos(resDf, filt):
  if filt == None:
    finalDf = resDf
  else:
    filtNa = resDf.loc[resDf['Var_Al_Relpos'] != 'NaN']
    filtNa['Var_Al_RelPos'] = filtNa['Var_Al_RelPos'].astype(float)
    finalDf = filtNa.loc[filtNa['Var_Al_RelPos'] >= filt]
  return(finalDf)
 
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
  try:
    resDf = isConsensus(df = resDf, freq = cons)
  except:
    print(f"Failed to find if variant is consensus, {current_directory}")
    sys.exit(1)
  try:
    resDf = adjFreq(df = resDf, relPos = rel)
  except:
    print(f"Failed to adjust frequencies {current_directory}")
    sys.exit(1)
  try:
    resDf = typeAllele(df = resDf)
  except:
    print(f"Failed to make corrected alleles and nucleotide change {current_directory}")
    sys.exit(1)
  ## annotations 
  # do we need these files as user input options?
  try:
    annot = pd.read_csv(f'{annot_dir}/py_annotations.csv')
    codons = pd.read_csv(f'{annot_dir}/programs/data/codons_table.csv', names = ['full_name', 'aa', 'codon', 'aa2', 'full_name2'])
  except:
    print(f"Failed to read annotation and/or codon tables {current_directory}")
    sys.exit(1)
  try:
    resDf = annotate(df=resDf, annot = annot, codons = codons)
  except:
    print(f"Failed annotating and translating mutations {current_directory}")
    sys.exit(1)
  try:
    resDf = ShPi(df = resDf)
  except:
    print(f"Failed to add Pi*Ln(Pi) for shannon index {current_directory}")
    sys.exit(1)
    
  resDf = filterResRelPos(resDf = resDf, filt = filt)
  return resDf

# execute runAll function
if __name__ == "__main__":
  resDf = runAll()
  writeResults(resDf = resDf, splitVar = splitVar, output = output)


    
  
