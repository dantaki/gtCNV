#!/usr/bin/env python
from src.preprocess import  preprocess
from src.feature_extraction import extract_feats
from src.genotype import genotype
import src.janitor as janitor
from src.vcf import annotate
import argparse
from argparse import RawTextHelpFormatter
import sys
from multiprocessing import Pool
import os
from glob import glob
from time import time
from datetime import timedelta 
class BAMlist:
        def __init__(self,fh):
                bam = {}
                gender = {}
                with open(fh) as f:
                        for l in f:
                                (id,bamfh,sex) = l.rstrip('\n').split('\t')
                                bam[id]=bamfh
                                gender[id]=sex
                self.bam=bam
                self.sex=gender
if __name__ == '__main__':
	init_time = int(time())
	splash ='        __________________   ____    __\n_______ __  /__  ____/__  | / /_ |  / /\n__  __ `/  __/  /    __   |/ /__ | / / \n_  /_/ // /_ / /___  _  /|  / __ |/ /  \n_\__, / \__/ \____/  /_/ |_/  _____/   \n/____/\nVersion 1.0\n\n\nAuthors: Danny Antaki <dantaki@ucsd.edu>, William Brandler\n'
        parser = argparse.ArgumentParser(description=splash,formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', help='Tab delimited input [ ID, BAM PATH, M/F ]',required=True)
	parser.add_argument('-b','--bed', help='Tab delimited BED file of CNVs [ CHROM, START, END, TYPE ]',type=str)
	parser.add_argument('-c','--cpu', help='Parallelize sample-wise. 1 per cpu',required=False,default=1,type=int)
        parser.add_argument('-g',  help='Reference genome build [ hg19, hg38 ]',required=False,default='hg19',type=str)
        parser.add_argument('-s','--seed', help='Preprocessing: integer seed for genome shuffling',required=False,default=42,type=int)
	parser.add_argument('-o','--out', help='outfile path',required=False,default="gtCNV_genotypes.vcf",type=str)
        parser.add_argument('--pre', help='Preprocessing output directory',required=False,default=None)
	parser.add_argument('--feats', help='Feature output directory',required=False,default=None)
	args = parser.parse_args()
        bam = args.i
        bed = args.bed
	cores = args.cpu
        gen = args.g
	ofh = args.out
        seed = args.seed
	predir= args.pre
	featsdir = args.feats
	preprocess_files={}
	feats_files={}
	janitor.check_in(bam)
	janitor.check_bed(bed)
	bam_dict = BAMlist(bam).bam
	gender_dict = BAMlist(bam).sex
	gens = ['hg19','hg38']
	if gen not in gens: 
		print "ERROR --genome must be either hg19 or hg38"
		sys.exit()
	"""----------------PREPROCESSING----------------"""
	if predir == None:
		if cores > 1 :
        		pool = Pool(processes=cores)
               		for bam_id in bam_dict:
				preofh = bam_id+'_gtCNV_preprocessing.txt'
              			preprocess_files[bam_id]=os.getcwd()+'/gtCNV_preprocessing/'+preofh
				pool.apply_async(preprocess, args=(bam_id,bam_dict[bam_id],preofh,gen,seed) )
                	pool.close()
                	pool.join()
		else:
        		for bam_id in bam_dict:
				preofh = bam_id+'_gtCNV_preprocessing.txt'
                        	preprocess_files[bam_id]=os.getcwd()+'/gtCNV_preprocessing/'+preofh
				preprocess(bam_id,bam_dict[bam_id],preofh,gen,seed)	
	else: 
		if not predir.endswith('/'): predir = predir+'/'
                if not os.path.isdir(predir):
                        print "ERROR "+predir+" PATH NOT FOUND"
                        sys.exit()
                for fh in glob(predir+'*gtCNV_preprocessing.txt'):
                        f = open(fh)
                        if sum(1 for l in open(fh)) <= 1: continue
                        else:
                                preids=[]
                                head=f.next().rstrip('\n')
                                if head == '\t'.join(("id","chr","cov","read_length_median","insert_size_median","insert_size_MAD","chr_bp_parsed")):
                                        for l in f: preids.append(l.rstrip('\n').split('\t').pop(0))
                        f.close()
                        for iid in set(preids): 
				if gender_dict.get(iid) != None: preprocess_files[iid]=fh
	print '-------- PREPROCESSING COMPLETE --------\n    time elapsed: '+str(timedelta(seconds=int(time())-init_time))+'\n----------------------------------------\n'
	""""----------------FEATURE EXTRACTION----------------"""
	if bed == None: 
		print "ERROR genotyping requires a BEDfile"
                sys.exit()
        if featsdir == None:
		if cores > 1 :
        		pool = Pool(processes=cores)
      			for bam_id in bam_dict:
				if preprocess_files.get(bam_id) == None:
					print 'ERROR preprocessing iid does not match inlist iid'
					continue
				prefh = preprocess_files[bam_id]
				gtofh = bam_id+'_gtCNV_features.txt'
				feats_files[bam_id]=os.getcwd()+'/gtCNV_genotypes/'+gtofh
				pool.apply_async(extract_feats, args=(bam_id,bam_dict[bam_id],bed,prefh,gender_dict[bam_id],gtofh,gen) )
                		pool.close()
                		pool.join()
		else:
			for bam_id in bam_dict:
				if preprocess_files.get(bam_id) == None:
                                	print 'ERROR preprocessing iid does not match inlist iid'
                                	continue
                        	prefh = preprocess_files[bam_id]
                        	gtofh = bam_id+'_gtCNV_features.txt'
                        	feats_files[bam_id]=os.getcwd()+'/gtCNV_genotypes/'+gtofh
				extract_feats(bam_id,bam_dict[bam_id],bed,prefh,gender_dict[bam_id],gtofh,gen)
	else: 
		if not featsdir.endswith('/'): featsdir=featsdir+'/'
		if not os.path.isdir(featsdir):
			print 'ERROR '+featsdir+' PATH NOT FOUND'
			sys.exit()
		for fh in glob(featsdir+'*gtCNV_features.txt'):
			f = open(fh)
                        if sum(1 for l in open(fh)) <= 1: continue
			else:
				featsid=[]
				head=f.next().rstrip('\n')
				if head == '\t'.join(('chr','start','end','type','size','id','coverage','discordant_ratio','split_ratio')):
					for l in f: featsid.append(l.rstrip('\n').split('\t').pop(5))
			f.close()
			for iid in set(featsid): 
				if gender_dict.get(iid) != None: feats_files[iid]=fh
	print '-------- FEATURE EXTRACTION COMPLETE --------\n    time elapsed: '+str(timedelta(seconds=int(time())-init_time))+'\n---------------------------------------------\n'
	"""----------------GENOTYPING----------------"""
	feats=[]
	for iid in feats_files:
		with open(feats_files[iid]) as f:
			for l in f: feats.append(tuple(l.rstrip('\n').split('\t')))
	genos,REF,NON,GQ,HEMI = genotype(feats,gender_dict,gen)
	annotate(genos,gen,REF,NON,GQ,ofh,gender_dict,HEMI)
	print '-------- GENOTYPING COMPLETE --------\n    time elapsed: '+str(timedelta(seconds=int(time())-init_time))+'\n-------------------------------------\n'
