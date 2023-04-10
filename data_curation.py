import pandas as pd
import glob
import re
import requests as r
from time import sleep
from Bio import SeqIO
from io import StringIO

path = "/Uniprot Data/"
filenames = glob.glob(path + "/*.xlsx")
print('No. of files:', len(filenames))

df_list = []
for file in filenames:
    file_path = file
    pattern = '[\w-]+?(?=\.)'
    a = re.search(pattern, file_path)
    x = str(a.group())
    print("Reading file = ",file)
    temp_df = pd.read_excel(file, sheet_name=x, header=None)
    temp_df.drop(temp_df.iloc[:, 1:2], inplace=True, axis=1)
    temp_df.drop(temp_df.iloc[:, 3:9], inplace=True, axis=1)
    temp_df.set_axis(['EntryID', 'PTM', 'Residue No.'], axis=1, inplace=True)
    temp_df['Modification'] = x
    df_list.append(temp_df)
    df = pd.concat(df_list, axis=0)
    df.reset_index(drop=True, inplace=True)
    df

df = pd.concat(df_list, axis=0)
df.reset_index(drop=True, inplace=True)

sequences = []
# c = 0
for cID in df['EntryID']:
  baseUrl="http://www.uniprot.org/uniprot/"
  currentUrl = baseUrl+str(cID)+".fasta"
  response = r.post(currentUrl)
  cData=''.join(response.text)
  Seq=StringIO(cData)
  pSeq=list(SeqIO.parse(Seq,'fasta'))
  seq = ''
  if(len(pSeq)!=0):
    seq = str(pSeq[0].seq)
    sequences.append(seq)
  # c+=1
  # print(c)

  sleep(0.03)

df['Sequence'] = sequences

# window size generation
mod_s = []
w = 7
for i in range(len(df)):
  x = df['Sequence'][i]
  y = df['Residue No.'][i]
  r = int(y+w) if ((y+w)<len(x)) else (len(x)-1)
  l = int(y-w-1) if ((y-w-1)>0) else 0
  # print((l, r))
  temp =  x[l:r]
  mod_s.append(temp)
len(mod_s)

df['Mod Sequence'] = mod_s

df = df[['EntryID',	'Sequence', 'Mod Sequence', 'PTM',	'Residue No.','Modification']]

df.to_csv('dataset.csv', index=False)