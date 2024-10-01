#!/usr/bin/env python
from Bio import SeqIO
import sys
from sys import argv
import re

#./add-read-lengths.py <samfile> <reads-calls2.txt>

samHandle=open('%s'%argv[1], 'r')
callHandle=open('%s'%argv[2],'r')
updateHandle=open('%s_added-lengths.txt'%argv[2],'w')
x=0
m=0
readDict={}
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
		if unmapped == 0 and secondary ==0 and fail ==0 and duplicate==0:
			seq=line.split()[9]
			if len(seq)>3000:
				m+=1
				readname=line.split()[0]
				mapcontig=line.split()[2]
				mappos=int(line.split()[3])
				tuple=[mappos, (mappos+len(seq))]
				#print('for read %s length is %s, %s, position %s, reversed %s'%(readname,len(seq),mapcontig,mappos, reversed))
				if readname in readDict.keys():
					oldtuple=readDict[readname]
					oldlen=(oldtuple[1]-oldtuple[0])
					if oldlen<len(seq):
						readDict[readname]=tuple #update the dictionary
						#print('duplicate entry for read %s, old len was %s new len %s, updated'%(readname, oldlen,len(seq)))
				else:
					readDict[readname]=tuple
		if x%1000==0:
			print('processed %s thousand reads from sam file'%(x/1000))
oldheader=next(callHandle)
hfirst=oldheader.split('\t')[:9]
updateHandle.write('\t'.join(hfirst))
updateHandle.write('\tmapped-segment-length\tcall-list\n')
j=0
for line in callHandle:
	j+=1
	first=line.split('\t')[:9]
	call_list=line.split('\t')[9]
	readname=line.split()[0]
	if readname in readDict.keys():
		tuple=readDict[readname]
		length=(tuple[1]-tuple[0])
		#print('for read %s tuple is %s and length is %s'%(readname,tuple,length))
		updateHandle.write('\t'.join(first))
		updateHandle.write('\t%s\t%s'%(length,call_list))
	if j%1000==0:
		print('processed %s thousand reads from call list file'%(j/1000))
samHandle.close()
callHandle.close()
updateHandle.close()
