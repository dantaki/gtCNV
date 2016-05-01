#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np
import glob
from sklearn.svm import SVC
import sys
from sklearn.externals import joblib
try:
   import cPickle as pickle
except:
   import pickle
def del_autosome_tab(df,preds,lik):
        df['copy_number']=preds
        df['likHOM'] = -10.0*np.log10(1.0-lik[:,0])
        df['likHET'] = -10.0*np.log10(1.0-lik[:,1])
        df['likREF'] = lik[:,2]
        df['NONREF'] = 1 - df['likREF']
        df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
        df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likHOM','likHET','PHRED_REF','PHRED_NONREF']]
	return df
def del_sexchr_tab (df,preds,lik):
        df['copy_number']=preds
        df['likVAR'] = -10.0*np.log10(1.0-lik[:,0])
        df['likREF'] = -10.0*np.log10(1.0-lik[:,1])
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likVAR','likVAR','likREF','likVAR']]
        return df
def dup_autosome_tab(df,preds,lik):
	df['copy_number']=preds
        df['likCN4'] = -10.0*np.log10(1.0-lik[:,0])
        df['likCN3'] = -10.0*np.log10(1.0-lik[:,1])
        df['likREF'] = lik[:,2]
	df['NONREF'] = 1 - df['likREF']
        df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
        df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
	df['copy_number'] = df['copy_number'].replace(0,4)
	df['copy_number'] = df['copy_number'].replace(1,3)
	df['copy_number'] = df['copy_number'].replace(2,2)	
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likCN4','likCN3','PHRED_REF','PHRED_NONREF']]
        return df
def dup_sexchr_tab(df,preds,lik):
	df['copy_number']=preds
        df['likVAR'] = -10.0*np.log10(1.0-lik[:,0])
        df['likREF'] = -10.0*np.log10(1.0-lik[:,1])
	df['copy_number'] = df['copy_number'].replace(0,2)
        df['copy_number'] = df['copy_number'].replace(1,1)
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likVAR','likVAR','likREF','likVAR']]
	return df
def merge3d (df):
        return np.vstack(np.asarray( (df['covr'],df['dpe'],df['sr']), order = 'C', dtype='float' )).T
def prep3d(df):
	df[['size']]=df[['size']].astype(int)
	df[['covr','dpe','sr']]=df[['covr','dpe','sr']].astype(float)
        df = df[df['covr'] < 5]
	X = merge3d(df)
	return df,X
def prep4d(df):
	df[['size']]=df[['size']].astype(int)
	df[['covr','dpe','sr']]=df[['covr','dpe','sr']].astype(float)
	df = df[df['covr'] < 5]
	big = df[df['size'] > 1000]
	sma = df[df['size'] <= 1000]
	X1 = merge3d(big)
	X2 = merge3d(sma)
	return big,sma,X1,X2
def autosome_del_svm(df):
	(bigdf,smadf,BIG,SMA) = prep4d(df)
	clf=''
	final=pd.DataFrame()
	if len(BIG) > 0:
		with open('resources/training_sets/1000genomes_gold_standard_gtCNV_svm_model_autosome_deletions_large_balanced_{}_{}.pkl'.format(0.01,10.0),'rb') as f: clf=pickle.load(f)
		preds = clf.predict(BIG)
		lik = clf.predict_proba(BIG)
		final=final.append(del_autosome_tab(bigdf,preds,lik))
	if len(SMA) > 0:
		clf=''
		with open('resources/training_sets/1000genomes_gold_standard_gtCNV_svm_model_autosome_deletions_small_balanced_{}_{}.pkl'.format(1000.0,1.0),'rb') as f: clf=pickle.load(f)
        	preds = clf.predict(SMA)
        	lik = clf.predict_proba(SMA)
        	final=final.append(del_autosome_tab(smadf,preds,lik))
	return final
def sexchr_del_svm(df):
	(bigdf,smadf,BIG,SMA) = prep4d(df)
        clf=''
        final=pd.DataFrame()
        if len(BIG) > 0:
		with open('resources/training_sets/1000genomes_gold_standard_gtCNV_svm_model_sexchr_deletions_large_weighted_{}_{}.pkl'.format(1.0,1.0),'rb') as f: clf=pickle.load(f)
        	preds = clf.predict(BIG)
        	lik = clf.predict_proba(BIG)
        	final=final.append(del_sexchr_tab(bigdf,preds,lik))
        if len(SMA) > 0:
		clf=''
	        with open('resources/training_sets/1000genomes_gold_standard_gtCNV_svm_model_sexchr_deletions_large_weighted_{}_{}.pkl'.format(1.0,1.0),'rb') as f: clf=pickle.load(f)
	        preds = clf.predict(SMA)
	        lik = clf.predict_proba(SMA)
	        final=final.append(del_sexchr_tab(smadf,preds,lik))
        return final
def autosome_dup_svm(df):
	(df,X) = prep3d(df)
        clf=''
	with open('resources/training_sets/1000genomes_gold_standard_gtCNV_svm_model_autosome_duplications_weighted_{}_{}.pkl'.format(1.0,10.0),'rb') as f: clf=pickle.load(f)
	preds = clf.predict(X)
        lik = clf.predict_proba(X)
        return dup_autosome_tab(df,preds,lik)
def sexchr_dup_svm(df):
	(df,X) = prep3d(df)
        clf=''
	with open('resources/training_sets/1000genomes_gold_standard_gtCNV_svm_model_sexchr_duplications_weighted_{}_{}.pkl'.format(1000.0,0.01),'rb') as f: clf=pickle.load(f)
	preds = clf.predict(X)
        lik = clf.predict_proba(X)
        return dup_sexchr_tab(df,preds,lik)	
