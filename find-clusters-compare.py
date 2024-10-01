#!/usr/bin/env python
from Bio import SeqIO
import sys
from sys import argv

#./find-clusters-compare.py <calls file 1> <integer 1 for depth> <calls file 2> <integer 2 for depth> 
x=int(argv[2])
y=int(argv[4])
countsDict={}
countsDict2={}
inHandle=open('%s'%argv[1],'r')
inHandle2=open('%s'%argv[3],'r')
next(inHandle)
next(inHandle2)
for line in inHandle:
	contig=line.split()[1]
	recomb=int(line.split()[5])
	position=int(line.split()[2])
	length=int(line.split()[9])
	stop=position+length
	tuple=(position,stop)
	if recomb==1:
		if contig in countsDict.keys(): 	
			previouslist=countsDict[contig]
			previouslist.append(tuple)
			countsDict[contig]=previouslist
			
		else:
			newlist=[]			
			newlist.append(tuple)
			countsDict[contig]=newlist

for line in inHandle2:
	contig=line.split()[1]
	recomb=int(line.split()[5])  
	position=int(line.split()[2])
	length=int(line.split()[9])
	stop=position+length 
	tuple=(position,stop)
	if recomb==1:
		if contig in countsDict2.keys():
			previouslist=countsDict2[contig]
			previouslist.append(tuple)
			countsDict2[contig]=previouslist
		else:
			newlist=[]
			newlist.append(tuple)
			countsDict2[contig]=newlist

runningtotal=0
contigs=countsDict.keys()
contigclusterlist1=[]

for c in contigs:
	tuplelist=countsDict[c]
	clusterCounter=1
	clusterstart=0
	newcluster=0
	totalclusters=0
	oldstart=500000000
	oldstop=0
	startlist=[]
	stoplist=[]
	#print('for %s list is %s'%(c,list))
	for tuple in tuplelist:
		newstart=int(tuple[0])
		newstop=int(tuple[1])
		if oldstart<newstart<oldstop:#overlapping
			clusterCounter+=1
			clusterflag=1
			startlist.append(newstart)
			stoplist.append(newstop)
		else: #newcluster starting?
			if clusterCounter>=x:
				clusterflag=0
				startlist.sort()
				cluststart=startlist[0]
				stoplist.sort(reverse=True)
				cluststop=stoplist[0]
				startlist.sort(reverse=True)
				maxstart=startlist[0]
				stoplist.sort()
				minstop=stoplist[0]
				totalclusters+=1
				runningtotal+=totalclusters
				if maxstart<minstop:
					distance=minstop-maxstart
					interval=(c,maxstart,minstop)
				else:#spread out
					distance=cluststop-cluststart
					interval=(c,cluststart,cluststop)
				contigclusterlist1.append(interval)

			startlist=[]
			stoplist=[]
			startlist.append(newstart)
			stoplist.append(newstop)
			clusterCounter=1
		oldstart=newstart
		oldstop=newstop
	if clusterflag==1: #unwritten data, rescue
		if clusterCounter>=x:
			startlist.sort()
			stoplist.sort(reverse=True)
			cluststart=startlist[0]
			cluststop=stoplist[0]
			startlist.sort(reverse=True)
			maxstart=startlist[0]
			stoplist.sort()    
			minstop=stoplist[0]
			if maxstart<minstop:
				distance=minstop-maxstart
				interval=(c,maxstart,minstop)
			else:
				interval=(c,cluststart,cluststop)#create tuple with minimal overlapping segment
				distance=cluststop-cluststart
			contigclusterlist1.append(interval)
			totalclusters+=1
			runningtotal+=totalclusters
		startlist=[]
		stoplist=[]
		clusterCounter=0
print('found %s clusters in file %s'%(len(contigclusterlist1),argv[1]))
runningtotal=0
contigs=countsDict2.keys()
contigclusterlist2=[]

for c in contigs:
	tuplelist=countsDict2[c]
	clusterCounter=1
	clusterstart=0
	newcluster=0
	totalclusters=0
	oldstart=500000000
	oldstop=0
	startlist=[]
	stoplist=[]
	for tuple in tuplelist:
		newstart=int(tuple[0])
		newstop=int(tuple[1])
		if oldstart<newstart<oldstop:#overlapping
			clusterCounter+=1
			clusterflag=1
			startlist.append(newstart)
			stoplist.append(newstop)
		else: #newcluster starting?
			if clusterCounter>=y:
				clusterflag=0
				startlist.sort()
				cluststart=startlist[0]
				stoplist.sort(reverse=True)
				cluststop=stoplist[0]
				startlist.sort(reverse=True)
				maxstart=startlist[0]
				stoplist.sort()
				minstop=stoplist[0]
				totalclusters+=1
				runningtotal+=totalclusters
				if maxstart<minstop:
					distance=minstop-maxstart
					interval=(c,maxstart,minstop)
				else:#spread out
					distance=cluststop-cluststart
					interval=(c,cluststart,cluststop)
				contigclusterlist2.append(interval)

			startlist=[]
			stoplist=[]
			startlist.append(newstart)
			stoplist.append(newstop)
			clusterCounter=1
		oldstart=newstart
		oldstop=newstop
	if clusterflag==1: #unwritten data, rescue
		if clusterCounter>=y:
			startlist.sort()
			stoplist.sort(reverse=True)
			cluststart=startlist[0]
			cluststop=stoplist[0]
			startlist.sort(reverse=True)
			maxstart=startlist[0]
			stoplist.sort()    
			minstop=stoplist[0]
			if maxstart<minstop:
				distance=minstop-maxstart
				interval=(c,maxstart,minstop)
			else:
				interval=(c,cluststart,cluststop)#create tuple with minimal overlapping segment
				distance=cluststop-cluststart
			contigclusterlist2.append(interval)
			totalclusters+=1
			runningtotal+=totalclusters
		startlist=[]
		stoplist=[]
		clusterCounter=0
print('found %s clusters in file %s'%(len(contigclusterlist2),argv[3]))
list1=[]
list2=[]
shared=[]
overlap=0
unique1=[]
unique2=[]
sharedDict={}
for cluster1 in contigclusterlist1:
	contig1=cluster1[0]
	start1=int(cluster1[1])
	stop1=int(cluster1[2])
	for cluster2 in contigclusterlist2:
		contig2=cluster2[0]
		start2=int(cluster2[1])
		stop2=int(cluster2[2])
		if contig1==contig2:
			if start1<start2<stop1 or start1<stop2<stop1 or start1<start2<stop2<stop1 or start2<start1<stop1<stop2:#overlapping
				overlap+=1
				combo=(cluster1,cluster2)
				shared.append(combo)
				sharedDict[cluster1]=1
				sharedDict[cluster2]=1
u1=0
u2=0
for cluster in contigclusterlist1:
	if cluster in sharedDict.keys():
		pass
	else:
		unique1.append(cluster)
		u1+=1

for cluster in contigclusterlist2:
	if cluster in sharedDict.keys():
		pass
	else:
		unique2.append(cluster)
		u2+=1
print('found %s unique clusters in file 1, %s'%(u1,argv[1]))
print('found %s unique clusters in file 2, %s'%(u2,argv[3]))
overl=(len(contigclusterlist1)-u1)
print('found %s overlapping clusters.'%overl)

print('file 1 clusters:')
for f1 in unique1:
	c1=f1[0]
	start=f1[1]
	stop=f1[2]
	print("%s\t%s\t%s"%(c1,start,stop))

print('**********')
#print('shared clusters:')		
for f1 in shared:
	c1=f1[0][0]
	start=f1[0][1]
	stop=f1[0][2]
	c2=f1[1][0]
	start2=f1[1][1]
	stop2=f1[1][2]
#	print("%s\t%s\t%s\t%s\t%s"%(c1,start,stop,start2,stop2))
#for f2 in unique2:
#	c2=f2[0]
#	start=f2[1]
#	stop=f2[2]
#	if c2=='Contig3':
#		print("%s\t%s\t%s"%(c2,start,stop))
		
inHandle.close()
inHandle2.close()
