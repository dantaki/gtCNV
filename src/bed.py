import pybedtools as pbed
class Bed:
	"""formats CNV input file as BED"""
	def __init__(self,fh):
		cnv=[]
		with open(fh) as f:
			for l in f:
				r = l.rstrip('\n').split('\t')
				cnv.append((r[0],r[1],r[2],r[3]))
		self.bed=pbed.BedTool(cnv).sort()
		temp=[]
		for i in self.bed: temp.append(i)
		self.cnv=temp
def annot_toList(cnv):
	cnv_list = []
	cnv.sort()
	for i in cnv:
		(c,s,e,cl,tag) = i
		cnv_list.append((c,s,e,cl,tag))
	return cnv_list
def bedsort_list(cnv):
	bed=pbed.BedTool(cnv).sort()
	cnv=[]
	for i in bed:
		(c,s,e,cl) = i
		cnv.append((c,s,e,cl))
	return cnv
def isPAR(cnv,gen):
	if len(pbed.BedTool(cnv,from_string=True).intersect(pbed.BedTool('resources/par_'+gen+'.bed'),f=0.5)) > 0: return True
	else: return False
def returnPAR(cnv,gen):
	out=[]
	results = pbed.BedTool(cnv).intersect(pbed.BedTool('resources/par_'+gen+'.bed'),f=0.5,wa=True,u=True)
	if len(results) > 0:
		for (c,s,e,cl,sz,iid,covr,dpe,sr) in results: out.append((c,s,e,cl,sz,iid,covr,dpe,sr))
	return out
def removePAR(cnv,gen):
	out=[]
	results= pbed.BedTool(cnv).subtract(pbed.BedTool('resources/par_'+gen+'.bed'),f=0.5,A=True)
	if len(results) > 0:
		for (c,s,e,cl,sz,iid,covr,dpe,sr) in results.sort(): out.append((c,s,e,cl,sz,iid,covr,dpe,sr))
	return out
