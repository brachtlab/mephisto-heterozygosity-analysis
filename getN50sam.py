#!/usr/bin/env python
from Bio import SeqIO
import sys
from sys import argv
import re

#./find-recombination.py <samfile> 

samHandle=open('%s'%argv[1], 'r')
x=0
lengthslist=[]
for line in samHandle:
	if line[:1] != "@": #not header
		x+=1
		readflag=int(line.split()[1])
		binstring=bin(readflag)
		reversebinstring=binstring[::-1]
		bpos=reversebinstring.find("b")
		binstringclip=reversebinstring[:bpos]
		length=len(binstringclip)
		flag=binstringclip.ljust(11,'0')
		unmapped=int(flag[2])
		reversed=int(flag[4])
		secondary=int(flag[8])
		fail=int(flag[9])
		duplicate=int(flag[10])
		if unmapped == 0:
			lengthslist.append(len(line.split()[9]))
		if x%10000==0:
			print('processed %s thousand reads.'%(x/1000))
runninglength=0
lengthslist.sort(reverse=True)
m=0
target=(sum(lengthslist)/2)
for length in lengthslist:
	m+=1
	runninglength+=length
	if runninglength>target:
		print('N50 is %s'%length)
		break
samHandle.close()
	

