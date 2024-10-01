#!/usr/bin/env python
from Bio import SeqIO
import sys
from sys import argv

#./add-context.py <variants_columns.txt> <genome.fasta>
genomeHandle=open('%s'%argv[2],'r')
snpHandle=open('%s'%argv[1],'r')
outHandle=open('%s_contextFIXED.txt'%argv[1],'w')
genomeDict={}
for seq_record in SeqIO.parse(genomeHandle,'fasta'):
	contig=seq_record.id
	seq=seq_record.seq
	genomeDict[contig]=seq

outHandle.write('%s\tprior7\tprior12\n'%(next(snpHandle).rstrip('\n')))
for line in snpHandle:
	ref=line.split()[7]
	alt=line.split()[8]
	contig=line.split()[0]
	pos=int(line.split()[1])
	prior7=genomeDict[contig][pos-8:pos-1]
	prior12=genomeDict[contig][pos-13:pos-1]
	if prior7!='':#ensure not empty
		if prior12!='':
			outHandle.write('%s\t%s\t%s\n'%(line.rstrip('\n'),prior7,prior12))

genomeHandle.close()
snpHandle.close()
outHandle.close()
