from random import sample
import pandas as pd
import os
# gtdb summary
df=pd.read_csv("data/gtdbtk.bac120.summary.tsv",sep='\t', converters={'user_genome': lambda x: str(x)})
# original metadata
metadat=pd.read_csv("data/metaDat.csv", converters={'uuid': lambda x: str(x)})
# gtdb NCBI dictionary
gtdbDat=pd.read_csv('C:/Users/abomb/OneDrive - Emory University/bacteria/data/gtdb/gtdbtkNcbiDict.tsv', sep='\t')

# get gtdbtk species names
def getGtdbSpec():
    phyDat=df[['classification']]
    splitPhyDat=phyDat['classification'].str.split(';', expand=True)
    splitPhyDatRen=splitPhyDat.rename(columns={splitPhyDat.columns[6]: 'Gtdb_taxa'})
    expTax=splitPhyDatRen[['Gtdb_taxa']]
    expTax['Gtdb_taxa']=expTax['Gtdb_taxa'].str.replace('s__', '')
    expTax['Gtdb_taxa']=expTax['Gtdb_taxa'].str.replace('Clostridium difficile', 'Clostridioides difficile')
    gtdbTax=pd.concat([df[['user_genome']], expTax], axis=1)
    gtdbTax = gtdbTax.rename(columns={'user_genome': 'uuid'})
    return gtdbTax

# get alternative taxa
def getAltTax():
    phyDat=df[['other_related_references(genome_id,species_name,radius,ANI,AF)']]
    splitPhyDat=phyDat['other_related_references(genome_id,species_name,radius,ANI,AF)'].str.split(',', expand=True)\
        .rename({ 1: 'Alt_taxa'}, axis=1)
    expTax=splitPhyDat[['Alt_taxa']]
    expTax['Alt_taxa']=expTax['Alt_taxa'].str.replace('s__', '')
    gtdbTax=pd.concat([df[['user_genome']], expTax], axis=1)
    gtdbTax = gtdbTax.rename(columns={'user_genome': 'uuid'})
    return gtdbTax

# get species names from original metadata
def getMetadatSpec():
    metadatSpec=metadat[['uuid', 'Organism']]
    metadatSpec=metadatSpec.rename(columns={'Organism':'Original_taxa'})
    metadatSpec['Original_taxa']=metadatSpec['Original_taxa'].str.replace('Clostridium difficile', 'Clostridioides difficile')
    return metadatSpec

# merge gtdbtk and original data
def mergeDat(df1, df2):
    mergedDat=pd.merge(df1, df2, left_on='uuid', right_on='uuid', how='left')
    return mergedDat

# merge main gtdb taxa and alternative taxa with metadata 
mainMerge=mergeDat(df1=getGtdbSpec(), df2=getMetadatSpec())
altMerge=mergeDat(df1=getAltTax(), df2=getMetadatSpec())

# get samples that match metadata
def findMatched():
    df=mainMerge
    matched=df.loc[df['Gtdb_taxa']==df['Original_taxa']]
    return matched

# get samples that mismatch metadata
def findMismatch():
    df=mainMerge
    mismatched=df.loc[df['Gtdb_taxa']!=df['Original_taxa']]
    return mismatched

def checkDictionary(dfTab, colName):
    dfTab['Ncbi_match']=''
    for i in range(0, len(dfTab.index)):
        # extract GTDB taxa name
        taxa=dfTab[colName][i]
        # match GTDB taxa name from results tables with dictionary
        dfMatch=gtdbDat[gtdbDat['gtdb_taxonomy'].str.contains(taxa)].reset_index()
        # if original name of extracted results column corresponds to dictionary GTDB to NCBI name write NCBI name
        if dfTab['Original_taxa'][i]==dfMatch['ncbi_organism_name'][0]:
            dfTab['Ncbi_match'][i]=dfMatch['ncbi_organism_name'][0]
        else:
            dfTab['Ncbi_match'][i]=='No_match'
    return dfTab

# extract species that mismatched due to spelling
def extractDictMatch():
    df=checkDictionary(dfTab=findMismatch(), colName='Gtdb_taxa')
    matched=df.loc[df['Ncbi_match']!='No_match']
    return matched

# extract taxa that should keep GTDB names
def matchRest():
    df=checkDictionary(dfTab=findMismatch(), colName='Gtdb_taxa')
    dfTab=df.loc[df['Ncbi_match']=='No_match']
    for i in range(0, len(dfTab.index)):
        taxa=dfTab['Gtdb_taxa'][i]
        dfMatch=gtdbDat[gtdbDat['gtdb_taxonomy'].str.contains(taxa)].reset_index()
        dfTab['Ncbi_match'][i]=dfMatch['ncbi_organism_name'][0]
    return dfTab

# merge results
def mergeResults():
    matched=findMatched()[['uuid', 'Gtdb_taxa']].reset_index()
    matched=matched.rename(columns={'Gtdb_taxa': 'taxa'})
    dictMatch=extractDictMatch()[['uuid', 'Original_taxa']].reset_index()
    rest=matchRest()[['uuid', 'Ncbi_match']]
    rest=rest.rename(columns={'Ncbi_match': 'taxa'})
    dictMatch=dictMatch.rename(columns={'Original_taxa': 'taxa'})
    finalTab=pd.concat([matched,dictMatch, rest], ignore_index=True).drop(['index'], axis=1)
    return finalTab

matchedResults=mergeResults()

#change directory
os.chdir('./test_gtdbtk')

# write bucket files and samples per taxa
taxaList=pd.unique(matchedResults['taxa'].tolist())
def writeSamples():
    for i in taxaList:
        #write name of the bacteria
        bucket=i.lower().replace(' ', '-')
        with open('./taxa/%s' % bucket, "a") as text_file:
            text_file.write(i + "\n")
        myData=matchedResults.loc[matchedResults['taxa']==i]
        samplesList=myData['uuid'].tolist()
        for j in samplesList:
            with open('./buckets/%s' % bucket, "a") as text_file:
                text_file.write(j + "\n")
               
writeSamples()

# check if all of the GTDB samples matched
def checkSamplesMatch():
    if len(matchedResults.index)==len(df.index):
        print('All results matched')
    else:
        data=matchRest()
        samples=data.loc[data['Ncbi_match']=='No_match']
        print('unmatched samples')
        print(samples)
        pd.DataFrame.to_csv(samples, 'unmatchedSamples.tsv', index=False, sep='\t')

checkSamplesMatch()





