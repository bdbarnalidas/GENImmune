# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 11:52:03 2021

@author: Maria del Carmen Martos contreras
@email: mcar.martos@pmcr.eu
"""
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import re
import argparse
import numpy as np


parser = argparse.ArgumentParser(prog="physiochemical_calculator", description='Peptide physiochemical calculator: Input must be a .csv file that inclundes the following columns names: peptide. This tool will generate different columns Molecular weight, instability, charge, GRAVY, Half life in the same file.  \nExample: \n  python3 phyisiochemical_calculator.py -f /path/to/peptide_pool.csv ')

parser.add_argument('--file', '-f', type=str,  help='csv file input with the peptides sequences to to calculate the physichochemical properties.')
args = parser.parse_args()

df = pd.read_csv (args.file, encoding='latin-1',sep=',', header=[0,1])
#df=df.drop_duplicates()
half_life = {'aa': ['A', 'R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'], 'mammal': [4.4,1,1.4,1.1,1.2,0.8,1,30,3.5,20,5.5,1.3,30,1.1,20,1.9,7.2,2.8,2.8,100]}
hl = pd.DataFrame(data=half_life)
def half_life(seq):
    index=hl.index
    cond = seq[0] == hl["aa"]
    ind=index[cond]
    return(hl["mammal"].loc[ind].item())

pi=[]
inst=[]
charg=[]
gravy=[]
half=[]
mw=[]

chars=["+", "(", ")", " ", ",", "-", "J", "X", "B", "6", "U", "7", "1", "2", "3", "4", "5", "8", "9", "l", "O", "Z", "\\", "Ë", "d", "i"] 

for i in df.index:
    pep=df[("Epitope", "Description")].loc[i]
    #pep=pep.encode('utf-8')
    if not ([extension for extension in  chars if(extension in pep)] or pep.islower()):
    	print(pep)
    	X = ProteinAnalysis(df[("Epitope", "Description")].loc[i])
    	mw.append("%0.2f" % X.molecular_weight())
    	pi.append("%0.2f" % X.isoelectric_point())
    	inst.append("%0.2f" % X.instability_index())
    	charg.append("%0.2f" % int(X.charge_at_pH(5.5)))
    	gravy.append("%0.2f" % X.gravy())    
    	half.append(half_life(df[("Epitope", "Description")].loc[i]))
    else:
    	#print(pep + " is not a readable peptide sequence")
    	mw.append(np.nan)
    	pi.append(np.nan)
    	inst.append(np.nan)
    	charg.append(np.nan)
    	gravy.append(np.nan)
    	half.append(np.nan)

df[('physicochemical','pI')],df[('physicochemical','Molecular_weight')],df[('physicochemical','instability')],df[('physicochemical','charge')],df[('physicochemical','GRAVY')],df[('physicochemical','Half_life')] = [pi,mw,inst,charg,gravy,half]


filter_df=df.dropna()

df.to_csv(args.file, sep=",")
filter_df.to_csv(args.file+"filtered.csv", sep=",")
