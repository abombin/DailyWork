import pandas as pd
import pytest
import os
#from itertools import islice
#import subprocess

#lofreq 2.1.5

def test_ReadCount(df1, df2):
  pyReadC = pd.read_table(df1, sep  = "\t")
  shReadC = pd.read_table(df2, sep = "\t")
  pyReadSort = pyReadC.sort_values(by = ['POS'])
  shReadSort = shReadC.sort_values(by = ['POS'])
  pyReadSort = pyReadSort.astype(str)
  shReadSort= shReadSort.astype(str)
  #assert pyReadSort.equals(shReadSort)
  if pyReadSort.equals(shReadSort):
    print("bam-readcount pass")
  else:
    print("bam-readcount FAIL")
    for i in range(len(shReadSort.index)):
      #print(i)
      position = pyReadSort.loc[i, "POS"]
      pyList = pyReadSort[pyReadSort["POS"] == position].values.flatten().tolist()
      shList = shReadSort[shReadSort["POS"] == position].values.flatten().tolist()
      if pyList != shList:
        #print(pyList)
        #print(shList)
        print(i)
        print('____________')
    

test_ReadCount(df1 = "sample_pos-filter.tsv" , df2 = "varcalltest_denv_04122023-1_varcall/position-filters/PG-DENV-056D_pos-filter.txt")

###
def test_lofreqFinal(df1, df2):
  pyReadC = pd.read_table(df1, sep  = "\t")
  shReadC = pd.read_table(df2, sep = "\t")
  pyReadSort = pyReadC.sort_values(by = ['POSITION'])
  shReadSort = shReadC.sort_values(by = ['POSITION'])
  pyReadSort['Pos_Ref_Alt'] = pyReadSort['POSITION'].astype(str) + pyReadSort['REF-NT'] + pyReadSort['VAR-NT']
  shReadSort['Pos_Ref_Alt'] = shReadSort['POSITION'].astype(str) + shReadSort['REF-NT'] + shReadSort['VAR-NT']
  pyReadFilt = pyReadSort[pyReadSort['Pos_Ref_Alt'].isin(shReadSort['Pos_Ref_Alt'])].sort_values(by = ['POSITION']).reset_index().drop(['index'], axis=1)
  #assert pyReadSort.equals(shReadSort)
  if pyReadFilt.equals(shReadSort):
    print("lowfreq final table pass")
  else:
    print("Initial lowfreq final table test FAIL")
    for i in range(len(pyReadFilt.index)):
      position = pyReadFilt.loc[i, "POSITION"]
      pyList = pyReadFilt[pyReadFilt["POSITION"] == position].values.flatten().tolist()
      shList = shReadSort[shReadSort["POSITION"] == position].values.flatten().tolist()
      if pyList != shList:
        print(i)
        print(pyList)
        print(shList)
        print('____________')
    
test_lofreqFinal(df1 = "sample_lofreq-output.tsv", df2 = "varcalltest_denv_04122023-1_varcall/final-calls/PG-DENV-056D_lofreq-output.txt")
