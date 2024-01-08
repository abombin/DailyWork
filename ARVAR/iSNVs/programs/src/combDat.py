import pandas as pd
import os
import glob

filesList = os.listdir('process')

def combDat(filesList, minFreq, maxFreq, minRelPos, maxSB):
  combDat = pd.DataFrame()
  for i in filesList:
    inFile = f'process/{i}/filtered.csv'
    try:
      curDf = pd.read_csv(inFile)
      filtDf = curDf[(curDf['ALLELE-FREQUENCY'] >= minFreq) & (curDf['ALLELE-FREQUENCY'] <= maxFreq) & (curDf['STRAND-BIAS'] <= maxSB) & (curDf['Var_Al_RelPos'] >= minRelPos)].reset_index().drop(["index"], axis =1)
      filtDf["Sample"] = i
      combDat = pd.concat([combDat, filtDf], axis = 0)
    except Exception as e:
      print(e)
  combDat = combDat.reset_index().drop(["index"], axis =1) 
  return combDat 

combDat = combDat(filesList = filesList, minFreq = 0.02, maxFreq = 1, minRelPos = 0.4, maxSB = 53)


def addUnSnv(df):
  Samp_Pos_Ref_Alt= []
  for i in (range(len(df.index))):
    curList = [df.loc[i, 'Sample'], str(df.loc[i, 'POSITION']), df.loc[i, 'REF-NT'], df.loc[i, 'VAR-NT']]
    curSamp = "__".join([ df.loc[i, 'Sample'], str(df.loc[i, 'POSITION']), df.loc[i, 'REF-NT'], df.loc[i, 'VAR-NT']])
    Samp_Pos_Ref_Alt.append(curSamp)
  df["Samp_Pos_Ref_Alt"] = Samp_Pos_Ref_Alt
  return df

combDat = addUnSnv(df = combDat)

combDat.to_csv("combDatFilt.csv", index = False)
