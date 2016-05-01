#!/usr/env python
import os
import sys
def check_in(fh):
	genders = ['M','F']
	if not os.path.isfile(fh):
		print 'ERROR '+fh+' not found'
		sys.exit()
	else:
		with open(fh) as f:
			for l in f:
				r = l.rstrip('\n').split('\t')
				if len(r) != 3: 
					print 'ERROR '+fh+' must be tab delimited\nFormat:\nIID    BAM file path   Gender[ M | F ]'
					sys.exit()
				else:
					(iid,bamfh,gender) = r
					if not os.path.isfile(bamfh):
						print 'ERROR '+bamfh+' not found'
						sys.exit()
					if gender not in genders:
						print 'ERROR accepted genders include M or F\n'+l
						sys.exit()
def check_bed(fh):
	chroms = [ 'chr1','chr2','chr3','chr4','chr5','chr6',
                   'chr7','chr8','chr9','chr10','chr11','chr12',
                   'chr13','chr14','chr15','chr16','chr17','chr18',
                   'chr19','chr20','chr21','chr22','chrX','chrY'
                 ]
	if not os.path.isfile(fh):
		print 'ERROR '+fh+' not found'
		sys.exit()
	else: 
		with open(fh) as f:
			for l in f:
				r = l.rstrip('\n').split('\t')
				c=r[0]
				s=r[1]
				e=r[2]
				cl=r[3]
				if c not in chroms: 
					print 'ERROR '+c+' not accepted chromosome format\nAccepted chromosomes:\n'
					print [k for k in chroms]
					sys.exit()
				if int(e) <= int(s):
					print 'ERROR '+e+' cannot be less than or equal to '+s+'\n'+l
					sys.exit()
				if 'DEL' not in cl:
					if 'DUP' not in cl:
						print 'ERROR '+cl+' must have DEL or DUP in the CNV type description\n'
						sys.exit()
