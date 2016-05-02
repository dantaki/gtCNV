#!/usr/env python
import pandas as pd
import numpy as np
import os
from glob import glob
from sklearn import svm
from math import ceil 
from inheritance import inheritance
from sklearn.preprocessing import StandardScaler
from bed import returnPAR,removePAR
from svm import autosome_del_svm,autosome_dup_svm,sexchr_del_svm,sexchr_dup_svm
try:
   import cPickle as pickle
except:
   import pickle
def genotype(feats,sex,gen):
	REF={}
	NON={}
	GQ={}
	HEMI={}
	males = [k for k in sex if sex[k] == 'M']
	sex_chrom = ['chrX','chrY']
	dels = [ k for k in feats if 'DEL' in k[3]]
	dups = [ k for k in feats if 'DUP' in k[3]]
	autosome_dels = [k for k in dels if k[0] not in sex_chrom]
        autosome_dups = [k for k in dups if k[0] not in sex_chrom]
        sexchr_dels = [k for k in dels if k[0] in sex_chrom]
        sexchr_dups = [k for k in dups if k[0] in sex_chrom]
	del_par=returnPAR(sexchr_dels,gen)
        dup_par=returnPAR(sexchr_dups,gen)
	if len(del_par) > 0:
                for k in del_par: autosome_dels.append(k)
        if len(dup_par) > 0:
                for k in dup_par: autosome_dups.append(k)
        male_sexchr_del=[]
        male_sexchr_dup=[]
        for k in removePAR(sexchr_dels,gen):
                if k[5] not in males: autosome_dels.append(k)
                else: male_sexchr_del.append(k)
        for k in removePAR(sexchr_dups,gen):
                if k[5] not in males: autosome_dups.append(k)
                else: male_sexchr_dup.append(k)
	head = ['chr','start','end','type','size','id','covr','dpe','sr']
	pd.options.mode.chained_assignment = None
	genos=[]
	if len(autosome_dels) > 0:
		autosome_del_df = pd.DataFrame(autosome_dels)
		autosome_del_df.columns=head
		for x in autosome_del_svm(autosome_del_df).values:
			x[6] = format(float(x[6])*2,'.2f')
			GQ[(x[0],x[1],x[2],x[4],x[5])]= ','.join((format(float(x[-4]),'.2f'),format(float(x[-3]),'.2f'),format(float(x[-2]),'.2f')))
			if int(x[9]) == 2: 
				x[9]='0/0'
				if x[5] not in males and x[0] == 'chrY': x[9]='.'
				else: 
					if REF.get((x[0],x[1],x[2],x[4]))==None: REF[(x[0],x[1],x[2],x[4])]=[float(x[-2])]
					else: REF[(x[0],x[1],x[2],x[4])].append(float(x[-2]))
			else:
				if int(x[9])==1: x[9]='0/1'
				else: x[9]='1/1'
				if x[5] not in males and x[0] == 'chrY': x[9]='.'
				else:
					if NON.get((x[0],x[1],x[2],x[4]))==None: NON[(x[0],x[1],x[2],x[4])]=[float(x[-1])]
                                	else: NON[(x[0],x[1],x[2],x[4])].append(float(x[-1]))
     			genos.append(tuple(x))
	if len(autosome_dups) > 0:
		autosome_dup_df = pd.DataFrame(autosome_dups)
        	autosome_dup_df.columns=head
		for x in autosome_dup_svm(autosome_dup_df).values:
			x[6] = format(float(x[6])*2,'.2f')
			GQ[(x[0],x[1],x[2],x[4],x[5])]= ','.join((format(float(x[-4]),'.2f'),format(float(x[-3]),'.2f'),format(float(x[-2]),'.2f')))
			if int(x[9]) == 2:
				x[9]='0/0'
				if x[5] not in males and x[0] == 'chrY': x[9]='.'
                                else: 
					if REF.get((x[0],x[1],x[2],x[4]))==None: REF[(x[0],x[1],x[2],x[4])]=[float(x[-2])]
                                	else: REF[(x[0],x[1],x[2],x[4])].append(float(x[-2]))
                        else:
				if int(x[9])==3: x[9]='0/1'
				else: x[9]='1/1'
				if x[5] not in males and x[0] == 'chrY': x[9]='.'
                                else:
					if NON.get((x[0],x[1],x[2],x[4]))==None: NON[(x[0],x[1],x[2],x[4])]=[float(x[-1])]
                                	else: NON[(x[0],x[1],x[2],x[4])].append(float(x[-1]))
			genos.append(tuple(x))
	if len(male_sexchr_del) > 0:
		sexchr_del_df = pd.DataFrame(male_sexchr_del)
        	sexchr_del_df.columns=head
		for x in sexchr_del_svm(sexchr_del_df).values:
			x[6] = format(float(x[6]),'.2f')
			GQ[(x[0],x[1],x[2],x[4],x[5])]= ','.join((format(float(x[-1]),'.2f'),format(float(x[-2]),'.2f')))
			HEMI[(x[0],x[1],x[2],x[4])]=1	
			if int(x[9]) == 1:
                                x[9]='0'
				if REF.get((x[0],x[1],x[2],x[4]))==None: REF[(x[0],x[1],x[2],x[4])]=[float(x[-2])]
                                else: REF[(x[0],x[1],x[2],x[4])].append(float(x[-2]))
                        else:
				x[9]='1'
				if NON.get((x[0],x[1],x[2],x[4]))==None: NON[(x[0],x[1],x[2],x[4])]=[float(x[-1])]
                                else: NON[(x[0],x[1],x[2],x[4])].append(float(x[-1]))
        		genos.append(tuple(x))
	if len(male_sexchr_dup) > 0:	
		sexchr_dup_df = pd.DataFrame(male_sexchr_dup)
        	sexchr_dup_df.columns=head
        	for x in sexchr_dup_svm(sexchr_dup_df).values:
			x[6] = format(float(x[6])*2,'.2f')
			GQ[(x[0],x[1],x[2],x[4],x[5])]= ','.join((format(float(x[-1]),'.2f'),format(float(x[-2]),'.2f')))
			HEMI[(x[0],x[1],x[2],x[4])]=1	
			if int(x[9]) == 1:
                                x[9]='0'
				if x[5] not in males and x[0] == 'chrY': x[9]='.'
				if REF.get((x[0],x[1],x[2],x[4]))==None: REF[(x[0],x[1],x[2],x[4])]=[float(x[-2])]
                                else: REF[(x[0],x[1],x[2],x[4])].append(float(x[-2]))
                        else:
				x[9]='1'
				if x[5] not in males and x[0] == 'chrY': x[9]='.'
                                if NON.get((x[0],x[1],x[2],x[4]))==None: NON[(x[0],x[1],x[2],x[4])]=[float(x[-1])]
                                else: NON[(x[0],x[1],x[2],x[4])].append(float(x[-1]))
			genos.append(tuple(x))
	refmed={}
	nonmed={}
	for x in REF: refmed[x]=int(ceil(np.median(REF[x])))
	for x in NON: nonmed[x]=int(ceil(np.median(NON[x])))
	return genos,refmed,nonmed,GQ,HEMI
