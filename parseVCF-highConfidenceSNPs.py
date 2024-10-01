#!/usr/bin/env python
import sys
from sys import argv
#./parseVCF2.py <vcf file> 
vcfHandle=open('%s'%argv[1],'r')
outHandle=open('%s_high-conf-snps.txt'%argv[1],'w')

outHandle.write("contig\tpos\ttotal_depth\tdepth-ref\tdepth-alt\tfraction-ref\tfraction-alt\tseq-ref\tseq-alt\tqual\n")
x=0
for line in vcfHandle:
	if line[:1]!='#':
		contig=line.split()[0]
		position=int(line.split()[1])
		vcf_ref=line.split()[3]
		vcf_alt=line.split()[4]
		qual=float(line.split()[5])
		info=line.split()[7]
		ilist=info.split(';')
		dp='null'
		dp4='null'
		refdepth='null'
		altdepth='null'
		reffrac='null'
		altfrac='null'
		for i in ilist:
			if 'DP=' in i:
				dp=i.split('=')[1]
			if 'DP4' in i:
				dp4=i.split('=')[1]
				#print(dp4)
				nums=dp4.split(',')
				refdepth=int(nums[0])+int(nums[1])
				altdepth=int(nums[2])+int(nums[3])
				totaldepth=refdepth+altdepth
				reffrac=float(refdepth)/(float(refdepth)+float(altdepth))
				altfrac=float(altdepth)/(float(refdepth)+float(altdepth))
		if qual>70 and altdepth > 4 and refdepth >4 and totaldepth<50:
			x+=1
			outHandle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(contig, position, (refdepth+altdepth), refdepth, altdepth, reffrac, altfrac, vcf_ref,vcf_alt,qual))

print("written %s sites"%x)
vcfHandle.close()
outHandle.close()
