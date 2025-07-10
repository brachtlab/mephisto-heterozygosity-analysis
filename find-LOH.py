#!/usr/bin/env python
import sys
from sys import argv
#./find-LOH.py <vcf1_columns.txt> #note, .txt 
					#files produced by parseVCF-freq.py
 
var1Handle=open('%s'%argv[1],'r')
outHandle=open('%s_LOH.txt'%(argv[1]),'w')
lengthsHandle=open('mephisto_contig_lengths.txt','r')
line1=next(var1Handle) #get past header
next(lengthsHandle)#get past header
lengthsDict={}
for line in lengthsHandle:
	contig=line.split()[0]
	length=int(line.split()[2])
	lengthsDict[contig]=length

lengthsHandle.close()

x=0
j=0
oldpos=0
oldcontig='null'
windowcounter=0
LOHflag=0
LOHstart=0
LOHstop=0
outHandle.write('contig\tstart\tstop\n')
for line in var1Handle:
	x+=1
	contig=line.split()[0]
	newpos=int(line.split()[1])
	if contig==oldcontig:
		gap=newpos-oldpos
		if gap>100000:
			LOHstart=oldpos
			LOHstop=newpos
			print('found intra-contig LOH: %s\t%s\t%s\tsize %s'%(contig,LOHstart,LOHstop,LOHstop-LOHstart))
			outHandle.write('%s\t%s\t%s\n'%(contig,LOHstart,LOHstop))
	else:#is new
		gap=newpos
		if oldcontig!='null':
			oldgap=lengthsDict[oldcontig]-oldpos
		else:
			oldgap=1
		if gap>100000:
			LOHstart=1
			LOHstop=newpos
			print('found contig-start LOH: %s\t%s\t%s\tsize %s'%(contig,LOHstart,LOHstop,LOHstop-LOHstart))
			outHandle.write('%s\t%s\t%s\n'%(contig,LOHstart,LOHstop))
		if oldgap>100000:#ran off end of last contig, rescue LOH from last contig
			end=lengthsDict[oldcontig]
			LOHstart=oldpos
			LOHstop=end
			print('found contig-end LOH: %s\t%s\t%s\tsize %s'%(oldcontig, LOHstart, LOHstop,LOHstop-LOHstart))
			outHandle.write('%s\t%s\t%s\n'%(oldcontig,LOHstart,LOHstop))
			LOHflag=0 #set back to zero
	oldcontig=contig
	oldpos=newpos

var1Handle.close()
outHandle.close()	

