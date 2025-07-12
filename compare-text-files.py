#!/usr/bin/env python
import sys
from sys import argv
#./compare-text-files.py <vcf1_columns.txt> <vcf2_columns.txt> #note, .txt 
					#files produced by parseVCF-freq.py
 
var1Handle=open('%s'%argv[1],'r')
var2Handle=open('%s'%argv[2],'r')
outHandle=open('compared_%s_&_%s'%(argv[1],argv[2]),'w')

Dict1={}
x=0
next(var1Handle)#get past header
for line in var1Handle:
	x+=1
	contig=line.split()[0]
	pos=line.split()[1]
	ident=contig+'_'+pos
	Dict1[ident]=line
print("number of variants in %s: %s"%(argv[1], x))

next(var2Handle)
both=0
one=0
two=0
j=0
m=0
for line2 in var2Handle:
	j+=1
	contig=line2.split()[0]
	pos=line2.split()[1]
	ident=contig+'_'+pos
	if ident in Dict1.keys():
		m+=1
		line1=Dict1[ident]
		outHandle.write("%s\t%s"%(line1.rstrip('\n'),line2))
		onefracalt=float(line1.split()[6])
		twofracalt=float(line2.split()[6])
		if onefracalt==1:
			if twofracalt==1:
				both+=1
			else:
				one+=1
		else:
			if twofracalt==1:
				two+=1
print("number of variants in %s: %s"%(argv[2],j))
print("found %s matching variants between the files, and of these:"%m)
print("number of homozygous alternate sites only in %s: %s"%(argv[1],one))
print("number of homozygous alternate sites only in %s: %s"%(argv[2],two))
print("number of homozygous sites in both: %s"%both)
			
var1Handle.close()
var2Handle.close()
outHandle.close()	
 

