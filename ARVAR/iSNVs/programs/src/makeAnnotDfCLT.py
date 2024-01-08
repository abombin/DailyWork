import time
start_time = time.time()
import argparse
import os
import pandas as pd
import numpy as np
import sys


parser=argparse.ArgumentParser(description='Make annotation data frame from grnr bank file')

# add arguments
parser.add_argument(
  "-i",
  "--input_path",
  type=str,
  help="Path to the input file, required"
)

parser.add_argument(
  "-o",
  "--output_path",
  type=str,
  nargs="?",
  default="py_annotations.csv",
  help="path to output file"
)

# parse arguments
args = parser.parse_args()

inFile = args.input_path
outFile = args.output_path

# functions

# read gb file
def getReadCountRes(inFile):
  with open(inFile, "r") as f:
    lines=f.readlines()
  return(lines)

# find 5' and 3' UTRs
def findUtrInd(genFile):
  for i in range(len(genFile)):
    if "5'UTR" in genFile[i]:
      seqStart = i
    if "3'UTR" in genFile[i]:
      seqEnd = i +1
  return seqStart, seqEnd

# parse gb file into lists
def ParseFile(seqStart, seqEnd, genFile):
  results = []
  for i in range(seqStart, seqEnd):
    newLine = genFile[i].split(" ")
    if newLine[5] == '':
      strLine = ' '.join(newLine)
      editLine = ' '.join(strLine.split())
    else:
      strLine = ' '.join(newLine)
      editLine = ' '.join(strLine.split())
      editLine = "__" + editLine
    results.append(editLine)

  resStr = ' '.join(results)
  editList = resStr.split("__")
  editList.pop(0)
  return(editList)

# get feature title and range
def getHeader(curList):
    headerList = curList[0].split(" ")
    curTitle = headerList[0]
    if "join(" not in headerList[1]:
      startPos, endPos = headerList[1].split("..")
    else:
      posList = headerList[1].replace("join(", "").replace(")", "").replace("..", ",").split(",")
      startPos = posList[0]
      endPos = posList[(len(posList) - 1)]
    return [curTitle, startPos, endPos]
  
# get title, positions, and features from the parsed list
def extractInfo(parsedList):
  allResults = []
  for i in range(len(parsedList)):
    curList = parsedList[i].split(" /")
    resList = getHeader(curList = curList)
    # extract elements from the list
    for elem in curList:
      if "gene=" in elem:
        curGene = eval(elem.replace("gene=", ""))
      elif "product=" in elem:
        curProduct = eval(elem.replace("product=", ""))
      elif "protein_id=" in elem:
        curProtId = eval(elem.replace("protein_id=", ""))
      elif "note=" in elem:
        curNote = eval(elem.replace("note=", ""))
      # deal with missing data  
      if 'curGene' not in locals():
        curGene = "NaN"
      if 'curProduct' not in locals():
        curProduct = "NaN"
      if 'curProtId' not in locals():
        curProtId = "NaN"
      if "curNote" not in locals():
        curNote = "NaN"
    
    resList.extend([curGene, curProduct, curProtId, curNote])
    del(curGene, curProduct, curProtId, curNote)
    allResults.append(resList)
    parsedResult = pd.DataFrame( allResults, columns = ["Title", "Start_pos", "End_pos", "Gene", "Product", "Protein_id", "Note"])
    parsedResult["Range"] = parsedResult["Start_pos"] + "-" + parsedResult["End_pos"]
    parsedResult[["Start_pos", "End_pos"]] = parsedResult[["Start_pos", "End_pos"]].apply(pd.to_numeric)  
    parsedResult["Length"] = parsedResult["End_pos"] - parsedResult["Start_pos"] +1
  return parsedResult

# get indicies for start and end of fasta sequence
def findSeq(genFile):
  for i in range(len(genFile)):
    if "ORIGIN" in genFile[i]:
      seqStart = i + 1
    if "//" in genFile[i]:
      seqEnd = i
  return seqStart, seqEnd

# get fasta sequence
def getFasta(fastStart, fastEnd, genFile):
  results = ""
  for i in range(fastStart, fastEnd, ):
    newLine = genFile[i].rstrip()
    editLine = ' '.join( newLine.split()).upper()
    lineList =  editLine.split(" ")
    filtLine = ''.join(lineList[1:])
    results = results + filtLine
  return results

# split features table ranges to cover all nucleotides individually
def splitFeatures(parsedResult):
  combDf = pd.DataFrame()
  for i in range(len(parsedResult.index)):
    newDf = pd.DataFrame(index=pd.Series(range(parsedResult.loc[i, "Start_pos"],parsedResult.loc[i, "End_pos"] +1), dtype = 'float64'))
    newDf["NT"] = newDf.index
    newDf["Title"] = parsedResult.loc[i, "Title"]
    newDf["Gene"] = parsedResult.loc[i, "Gene"]
    newDf["Product"] = parsedResult.loc[i, "Product"]
    newDf["Protein_id"] = parsedResult.loc[i, "Protein_id"]
    newDf["Length"] = parsedResult.loc[i, "Length"]
    combDf = pd.concat([combDf, newDf], ignore_index=True)
  sortDf = combDf.sort_values(by=['NT']).reset_index().drop(['index'], axis = 1)
  return(sortDf)

# filter NaN from the list
def filterNA(x):
  newList = []
  for i in x:
    if i != 'NaN':
      newList.append(i)
  return newList

# get features per position that have the smallest sequence length supporting the next function
def getMinVar(positionList, header, nuclRef):
  for position in positionList:
      minRef = nuclRef[nuclRef["Length"] == position]
      combVar = list(np.unique(minRef[header]))
      if any(val != "NaN" for val in combVar):
        combVar = filterNA(combVar)
        combStr = ",".join(combVar)
        break
  if 'combStr' not in locals():
    combStr = "NaN"
  return combStr

# get fearues per nucleotide for a new annotation table
def annotateNucl(refDf, i):
  nuclRef = refDf[refDf["NT"] == i]
  if len(nuclRef.index) > 0:
    nuclRef = nuclRef.sort_values(by = ["Length"])
    headList = ["Title", "Gene", "Product", "Protein_id"]
    positionList = list(np.unique(nuclRef["Length"]))
    refList = []
    minList = []
    for header in headList:
      combVar = np.unique(nuclRef[header])
      combStr = ",".join(filterNA(combVar))
      if combStr == '':
        combStr = "NaN"
      refList.append(combStr)
    for header in headList:
      minVar = getMinVar(positionList, header, nuclRef)
      minList.append(minVar)
  else:
    refList = ["NCR"] *4
    minList = ["NCR"] *4
  return refList, minList

# run previous two functions across all NT positions and add features to annotation table
def annotFastaDf(refDf, fastaDf):
  nuclPositions = fastaDf["NT"].to_list()
  allFullRef = []
  allMinRef = []
  for i in nuclPositions:
    fullRef, minRef = annotateNucl(refDf, i)
    allFullRef.append(fullRef)
    allMinRef.append(minRef)
  
  fullRefDf = pd.DataFrame(allFullRef, columns = ["All_Titles", "All_Genes", "All_Products", "All_Protein_IDs"])
  minRefDf = pd.DataFrame(allMinRef, columns = ["Min_Titles", "Min_Genes", "Min_Products", "Min_Protein_IDs"])
  annotDf = pd.concat([fastaDf, minRefDf, fullRefDf], axis = "columns")
  annotDf["All_Products_Edit"] = reformPrevNext(my_list = annotDf["All_Products"].to_list(), n = "NCR")
  return(annotDf)

# find previous element after n in the list that is not n
def prevElem(my_list, i, n):
  for elem in range(i-1, -1, -1):
    prevElem = my_list[elem]
    if prevElem != n:
      newElem = prevElem
      break
  return newElem

# find next element after n in the list that is not n
def nextElem(my_list, i, n):
  for elem in range(i, len(my_list)):
    nextElem = my_list[elem]
    if nextElem != n:
      newElem = nextElem
      break
  return newElem

# reannotate NCRs based on previous and next elements
def reformPrevNext(my_list, n):
  resList = []
  for i in range(len(my_list)):
    curVal = my_list[i]
    if curVal == n:
      prevVal = prevElem(my_list, i, n)
      nextVal = nextElem(my_list, i, n)
      newVal = "_".join([str(prevVal), str(curVal), str(nextVal)])
    else:
      newVal = str(curVal)
    resList.append(newVal)
  return resList

# find position of the first codon after 5' UTR
def findStartCod(annotDf):
  for i in range(len(annotDf.index)):
    if "5'UTR" not in annotDf.loc[i, "All_Titles"]:
      targNuc = annotDf.loc[i, "NT"]
      break
  return targNuc

# find position of the last codon before 3' UTR
def findEndCod(annotDf):
  for i in range(len(annotDf.index)):
    if "3'UTR" in annotDf.loc[i, "All_Titles"]:
      targNuc = annotDf.loc[i-1, "NT"]
      break
  return targNuc

# add codons per nucleotide in annotation table
def addCodons(fastaSeq, startCod, endCod, annotDf ):
  utr5Seq = ["NCR"] * (startCod-1)
  utr3Seq = (len(fastaSeq) - endCod) * ["NCR"]
  subStr = fastaSeq[(startCod -1):endCod]
  # split string in a list of 3 elements
  substrList = [subStr[i:i+3] for i in range(0, len(subStr), 3)]
  cdsList = []
  for item in substrList:
    cdsList.extend([item]*3)
  seqList = utr5Seq + cdsList + utr3Seq
  len(seqList) == len(fastaSeq)
  annotDf["Ref_Codon"] = seqList
  return(annotDf)

# count codon number per feature represented by all products overlappin with nucleotide
def addCodNumbPerBase(annotDf):
  editDf = pd.DataFrame()
  productsList = list(pd.unique(annotDf["All_Products_Edit"]))
  for i in productsList:
    curProd = annotDf[annotDf['All_Products_Edit'] == i].reset_index().drop(['index'], axis =1)
    curProd["Codon#"] = "NaN"
    for j in range(len(curProd.index)):
      if (curProd.loc[j, "All_Products"] == "NaN") or (curProd.loc[j, "All_Products"] == "NCR"):
        curProd.loc[j, "Codon#"] = "NCR"
      else:
        if j < 3:
          curProd.loc[j, "Codon#"] = 1
        else:
          newCod = curProd.loc[j-3, "Codon#"] + 1
          curProd.loc[j, "Codon#"] = newCod
    editDf = pd.concat([editDf, curProd])
  editDfSort = editDf.sort_values(by = ["NT"]).reset_index().drop(["index"], axis =1)
  return(editDfSort)

# translate codons to amino acids
def transLateCod(annotDf):
  aaList = []
  codons = pd.read_csv('programs/data/codons_table.csv', names = ['full_name', 'aa', 'codon', 'aa2', 'full_name2'])
  for i in range(len(annotDf.index)):
    codon = annotDf.loc[i, "Ref_Codon"]
    if codon == "NCR" or codon == "NaN":
      newAA = "NCR"
    else:
      newAA = codons.loc[codons['codon'] == codon, 'aa'].iloc[0]
    aaList.append(newAA)
  annotDf["Ref_AA"] = aaList
  return(annotDf)

# save annotation file
def saveAnnot(annotDf, outFile):
  if len(annotDf.index) > 0:
    annotDf = annotDf.rename(columns = {'Min_Products' : 'Region/Gene'})
    annotDf.to_csv(path_or_buf = outFile, index = False)

# combine all functions to make annotation table
def runAll(inFile):
  try:
    genFile = getReadCountRes(inFile)
    seqStart, seqEnd = findUtrInd(genFile)
    parsedFile = ParseFile(seqStart, seqEnd, genFile)
  except:
    print("Failed to read file")
    sys.exit(1)
  try:
    parsedResult = extractInfo(parsedList = parsedFile)
  except:
    print("Parsing file failed")
    sys.exit(1)
  try:
    fastStart, fastEnd = findSeq(genFile)
    fastaSeq = getFasta(fastStart, fastEnd, genFile)
    fastaDf = pd.DataFrame(list(fastaSeq), columns =['Ref'])
    fastaDf.insert(loc = 0, column = "NT", value = fastaDf.index +1)
  except:
    print("Extraction of the fasta sequence failed")
    sys.exit(1)
  try:
    refDf = splitFeatures(parsedResult)
  except:
    print("Features failed to be recorded as data frame")
    sys.exit(1)
  try:
    annotDf = annotFastaDf(refDf, fastaDf)
  except:
    print("Annotations per nucleotide failed")
    sys.exit(1)
  try:
    startCod = findStartCod(annotDf)
    endCod = findEndCod(annotDf)
    annotDf = addCodons(fastaSeq, startCod, endCod, annotDf)
  except:
    print("Failed adding codons for each nucleotide")
    sys.exit(1)
  try:
    annotDf = addCodNumbPerBase(annotDf)
  except:
    print("Failed to add codon numbers")
    sys.exit(1)
  try:
    annotDf = transLateCod(annotDf)
  except:
    print("Failed to translate codons")
    sys.exit(1)
  return(annotDf)

# run the tool
if __name__ == "__main__":
  annotDf = runAll(inFile)
  saveAnnot(annotDf = annotDf , outFile = outFile)

print("--- %s minutes ---" % ((time.time() - start_time) /60) )
