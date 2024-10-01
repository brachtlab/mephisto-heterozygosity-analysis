#!/usr/bin/env python
import sys
from sys import argv
#./filter-text-files.py <vcf1_columns.txt> #note, .txt 
					#files produced by parseVCF-freq.py
 
var1Handle=open('%s'%argv[1],'r')
outHandle=open('%s_snps_only.txt'%(argv[1]),'w')

line1=next(var1Handle) #get past header
outHandle.write(line1)
x=0
j=0
for line in var1Handle:
	x+=1
	ref=line.split()[7]
	alt=line.split()[8]
	if len(ref) == 1 and len(alt)==1:
		outHandle.write(line)
		j+=1
print("of %s entries, wrote %s snp-only."%(x,j))
var1Handle.close()
outHandle.close()	
 

