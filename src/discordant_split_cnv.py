#!/usr/env python
import pysam
def is_discordant(read,windows,ci,matepos):
        """returns True if read is discordant"""
        ((s1,e1),(s2,e2)) = windows
        if abs(read.tlen) >= ci:
                """insert size is greater than 5MAD from median"""
                if (int(s1) <= read.pos+1 <= int(e1) and int(s2) <= int(matepos)+1 <= int(e2)) or (int(s2) <= read.pos+1 <= int(e2) and int(s1) <= int(matepos)+1 <= int(e1)): 
			return True
                else:
                        return False
def is_split(read,windows,c):
        """returns True if read is split"""
        ((s1,e1),(s2,e2)) = windows
        if read.is_secondary == True:
                second_align = read.get_tag("SA").split(',')
                if second_align[0] != c:
                        """return False if maps to other chromosome"""
                        return False
                else:
                        second_align[1] = int(second_align[1])
                        if ( ((int(s1) <= read.pos+1 <= int(e1)) and (int(s2) <= second_align[1] <= int(e2))) 
			     or ((int(s2) <= read.pos+1 <= int(e2)) and (int(s1) <= second_align[1] <= int(e1)))):
                                """secondary alignment must be near the opposite breakpoint of the primary alignment"""
                                return True
                        else:
                                return False
def discordant_split_cnv(flank_list,bam,size,ci,windows,chrFlag):
        """count discordant paired ends and split reads"""
        discordant_count=0
        split_count=0
        concordant_count=0
        if flank_list==None: return (0,0,0)
	else: 
		for (c,s,e) in flank_list:
                	if chrFlag == False : c = c.replace("chr","")
                	region = str(c+":"+s+"-"+e)
                	for read in bam.fetch(region=region,until_eof=True):
                        	"""skip noninformative reads
                        	- courtsey of svtyper - https://github.com/hall-lab/svtyper --"""
                        	if (read.is_qcfail
                        	    or read.is_unmapped
                        	    or read.mate_is_unmapped
                        	    or read.tid != read.rnext
				    or read.is_reverse == read.mate_is_reverse
                        	    or read.is_duplicate): continue
                        	else:
					mate_pos = read.pnext
					if is_discordant(read,windows,ci,mate_pos) == True: discordant_count+=1
					if is_split(read,windows,c)  == True: split_count+=1
                                	if (read.is_proper_pair
                                	    and read.mapping_quality >= 10
					    and not read.is_secondary
                                	    and abs(read.tlen) < ci):
                                	        concordant_count+=1
                                	else: continue
        	return(discordant_count,split_count,concordant_count)
