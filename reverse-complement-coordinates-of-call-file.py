#!/usr/bin/env python
from Bio import SeqIO
import sys
from sys import argv

#./reverse-complement-coordinates-of-call-file.py <variants_columns.txt> <genome.fasta>
genomeHandle=open('%s'%argv[2],'r')
callHandle=open('%s'%argv[1],'r')
outHandle=open('%s_rc.txt'%argv[1],'w')
genomeDict={}
for seq_record in SeqIO.parse(genomeHandle,'fasta'):
	contig=seq_record.id
	seq=seq_record.seq
	genomeDict[contig]=len(seq)
#outHandle.write('readname\tcontig\tmap-position\treference-flag\talternate-flag\trecombination-flag\treference-counts\talt-counts\terror-counts\tcall-list\n')
next(callHandle)
for line in callHandle:
	readname=line.split()[0]
	contig=line.split()[1]
	refflag=line.split()[3]
	altflag=line.split()[4]
	recombflag=line.split()[5]
	refcounts=line.split()[6]
	altcounts=line.split()[7]
	errcounts=line.split()[8]
	position=int(line.split()[9])
	calllist=line.split('\t')[10].rstrip('\n')
	correctposition=(genomeDict[contig]-position)
	outHandle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(readname,contig.rstrip('_rc'),correctposition,refflag,altflag,recombflag,refcounts,altcounts,errcounts,calllist))

genomeHandle.close()
callHandle.close()
outHandle.close()
