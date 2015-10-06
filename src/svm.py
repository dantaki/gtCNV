#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np
import glob
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
import sys
import pickle
from sklearn.externals import joblib
######################################
def del_autosome_tab (df,preds,lik):
        df['copy_number']=preds
        df['likHOM'] = lik[:,0]
        df['likHET'] = lik[:,1]
        df['likREF'] = lik[:,2]
        df['NONREF'] = 1 - df['likREF']
        df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
        df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likREF','NONREF','PHRED_REF','PHRED_NONREF']]
        return df
def del_sexchr_tab (df,preds,lik):
        df['copy_number']=preds
        df['likHOM'] = lik[:,0]
        df['likHET'] = "NA"
        df['likREF'] = lik[:,1]
        df['NONREF'] = 1 - df['likREF']
        df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
        df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likREF','NONREF','PHRED_REF','PHRED_NONREF']]
        return df
def dup_autosome_tab(df,preds,lik):
	df['copy_number']=preds
        df['likREF'] = lik[:,0]
        df['likGAIN'] = lik[:,1]
        df['NONREF'] = 1 - df['likREF']
        df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
        df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
	df['copy_number'] = df['copy_number'].replace('0','2')
	df['copy_number'] = df['copy_number'].replace('1','3')	
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likREF','NONREF','PHRED_REF','PHRED_NONREF']]
        return df
def dup_sexchr_tab(df,preds,lik):
	df['copy_number']=preds
        df['likREF'] = lik[:,0]
        df['likGAIN'] = lik[:,1]
        df['NONREF'] = 1 - df['likREF']
        df['PHRED_REF'] = -10.0 * np.log10(1.0-df['likREF'])
        df['PHRED_NONREF'] = -10.0 * np.log10(1.0-df['NONREF'])
        df['copy_number'] = df['copy_number'].replace('1','2')
        df['copy_number'] = df['copy_number'].replace('0','1')
        df= df[['chr','start','end','size','type','id','covr','dpe','sr','copy_number','likREF','NONREF','PHRED_REF','PHRED_NONREF']]
        return df
def merge3d (df):
        return np.vstack(np.asarray( (df['covr'],df['dpe'],df['sr']), order = 'C', dtype='float' )).T
def svm_prepare(df):
	df.columns = ["chr","start","end","size","type","id","covr","dpe","sr"]
        df = df.fillna(0)
        X = merge3d(df)
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
	return df,X
def autosome_del_svm(df):
	(df,X) = svm_prepare(df)
	clf = joblib.load('resources/training_sets/svm/1kg_svm_model_autosome_del.pkl')
	preds = clf.predict(X)
	lik = clf.predict_proba(X)
	return del_autosome_tab(df,preds,lik)
def sexchr_del_svm(df):
	(df,X) = svm_prepare(df)
	clf = joblib.load('resources/training_sets/svm/1kg_svm_model_sexchr_del.pkl')
        preds = clf.predict(X)
        lik = clf.predict_proba(X)
        return del_sexchr_tab(df,preds,lik)
def autosome_dup_svm(df):
	(df,X) = svm_prepare(df)
	clf = joblib.load('resources/training_sets/svm/1kg_svm_model_autosome_dup.pkl')
        preds = clf.predict(X)
        lik = clf.predict_proba(X)
        return dup_autosome_tab(df,preds,lik)
def sexchr_dup_svm(df):
	(df,X) = svm_prepare(df)
	clf = joblib.load('resources/training_sets/svm/1kg_svm_model_sexchr_dup.pkl')
        preds = clf.predict(X)
        lik = clf.predict_proba(X)
        return dup_sexchr_tab(df,preds,lik)	
