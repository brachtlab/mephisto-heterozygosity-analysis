#!/usr/bin/env python
from Bio import SeqIO
import sys
from sys import argv
import re

#./find-recombination.py <samfile> <variants-file.vcf>

samHandle=open('%s'%argv[1], 'r')
x=0
m=0
recombcounter=0
callHandle=open('%s_reads_calls2.txt'%argv[1],'w')
callHandle.write('readname\tcontig\tmap-position\treference-flag\talternate-flag\trecombination-flag\treference-counts\talt-counts\terror-counts\tmappos+readlen\tcall-list\n')
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
		if unmapped == 0 and secondary ==0 and fail ==0 and duplicate==0 and reversed==0:
			seq=line.split()[9]
			if len(seq)>3000:
				m+=1
				readname=line.split()[0]
				mapcontig=line.split()[2]
				mappos=int(line.split()[3])
				cigar=line.split()[5]
				foundcontig=0
				varHandle=open('%s'%argv[2],'r')
				next(varHandle)
				vlist=[]
				v=0
				refcounter=0
				altcounter=0
				errorcounter=0
				call_list=[]
				for l in varHandle:
					vcontig=l.split()[0]
					vpos=int(l.split()[1])
					if vcontig==mapcontig:
						foundcontig=1
						if mappos <=vpos<=(mappos+len(seq)):
							vlist.append(l)
							v+=1
					else:
						if foundcontig==1:
							break #we went past contig in var file				
				varHandle.close()
				splitcigar=cigar.split('M')
				#print(splitcigar)
				vc=0
				cigarcounter=0
				runningdel=0
				runningins=0
				pointer=0
				fs=0
				found=0
				if 'S' in splitcigar[0]:
					softclip=int(splitcigar[0].split('S')[0])
				else:
					softclip=0
				adjustvarpos=1
				for var in vlist:#convert to read-space
					vc+=1
					vacontig=var.split()[0]
					vapos=int(var.split()[1])-(mappos-softclip)
					if vapos>=len(seq):
						break
					varef=var.split()[7]
					vaalt=var.split()[8]
					genomeprevious=var.split()[10]
					longprevious=var.split()[11]
					if pointer >=(vapos+runningins-runningdel):
						readvar=seq[vapos+runningins-runningdel]
						readprevious=seq[(vapos+runningins-runningdel-8):(vapos+runningins-runningdel-1)]
						if genomeprevious==readprevious:
							#print('-->1variant: reference %s alternate %s, read %s,  genomeprevious %s, readprevious %s, read-coord-of-var is %s, pointer %s, running deletions is %s and running insertions is %s'%(varef,vaalt, readvar, genomeprevious, readprevious, (vapos+runningins-runningdel) ,pointer, runningdel, runningins))
							found+=1
							if readvar==varef:
								#print('reference')
								call_list.append('ref')
								refcounter+=1
							elif readvar==vaalt:
								#print('alt')
								call_list.append('alt')
								altcounter+=1
							else:#is error
								errorcounter+=1
								#call_list.append('error')
						else:
							plusone=seq[(vapos+runningins-runningdel-7):(vapos+runningins-runningdel)]
							plustwo=seq[(vapos+runningins-runningdel-6):(vapos+runningins-runningdel+1)]
							plusthree=seq[(vapos+runningins-runningdel-5):(vapos+runningins-runningdel+2)]
							minusone=seq[(vapos+runningins-runningdel-9):(vapos+runningins-runningdel-2)]
							minustwo=seq[(vapos+runningins-runningdel-10):(vapos+runningins-runningdel-3)]
							minusthree=seq[(vapos+runningins-runningdel-11):(vapos+runningins-runningdel-4)]
							if genomeprevious==plusone:
								#print('found it shifted +1')
								context=seq[(vapos+runningins-runningdel-7):(vapos+runningins-runningdel+1)]
							elif genomeprevious==plustwo:
								#print('found it shifted +2')
								context=seq[(vapos+runningins-runningdel-6):(vapos+runningins-runningdel+2)]
							elif genomeprevious==plusthree:
								#print('found it shifted +3')
								context=seq[(vapos+runningins-runningdel-5):(vapos+runningins-runningdel+3)]
							elif genomeprevious==minusone:
								#print('found it shifted -1')
								context=seq[(vapos+runningins-runningdel-9):(vapos+runningins-runningdel-1)]
							elif genomeprevious==minustwo:
								#print('found it shifted -2')
								context=seq[(vapos+runningins-runningdel-10):(vapos+runningins-runningdel-2)]
							elif genomeprevious==minusthree:
								#print('found it shifted -3')
								context=seq[(vapos+runningins-runningdel-11):(vapos+runningins-runningdel-3)]
							else:
								#print('could not be found within 6 bp.')
								match=seq.find(longprevious)
								if match != -1:
									context=seq[match:match+13]
									#print('longprevious is %s and context is %s and base is %s found at index %s'%(longprevious, context,context[::-1][0], match))
								else:
									readvar='null'
							if readvar !='null':
								readvar=context[::-1][0]#get last item
								#print('-->1variant: reference %s alternate %s, read %s,  genomeprevious %s, readprevious %s, read-coord-of-var is %s, pointer %s, running deletions is %s and running insertions is %s'%(varef,vaalt, readvar, genomeprevious, readprevious, (vapos+runningins-runningdel) ,pointer, runningdel, runningins))
								found+=1
								if readvar==varef:
									#print('reference')
									call_list.append('ref')
									refcounter+=1
								elif readvar==vaalt:
									#print('alt')
									call_list.append('alt')
									altcounter+=1
								else:
									errorcounter+=1	
									#call_list.append('error')							
					else: #need to catch up
						while pointer < (vapos+runningins-runningdel): #the goal here is to map the indels leading up the read position of the variant. Pointer tracks CIGAR progress.
							if  cigarcounter+1 >=len(splitcigar):
								break
							c=splitcigar[cigarcounter]
							#print('cigar %s is number %s'%(c,cigarcounter+1))
							#print('in while loop pointer is %s and read-va-pos is %s'%(pointer,vapos+runningins-runningdel))
							if 'S' in c:
								#print('fs is %s'%fs)
								if fs==0:#first soft-clip
									soft=int(c.split('S')[0])
									match=int(c.split('S')[1])
									pointer=match+soft
									fs=1
								elif fs==1: #end of read
									break
							elif 'D' in c:
								if c.split('D')[0].isnumeric() and c.split('D')[1].isnumeric():
									delete=int(c.split('D')[0])
									match=int(c.split('D')[1])
									pointer-=delete
									pointer+=match
									runningdel+=delete
							elif 'I' in c:
								if c.split('I')[0].isnumeric() and c.split('I')[1].isnumeric():
									insert=int(c.split('I')[0])
									match=int(c.split('I')[1])
									pointer+=insert
									pointer+=match
									runningins+=insert
							cigarcounter+=1
						#print('readname %s seqlen %s index %s runningins %s runningdel %s readpos %s'%(readname,len(seq),vapos+runningins-runningdel,runningins,runningdel,vapos))
						#print('readlen %s vapos %s runningins %s runningdel %s'%(len(seq),vapos,runningins,runningdel))
						if (vapos+runningins-runningdel)>=len(seq) or (vapos+runningins-runningdel)<=10:
							break
						readvar=seq[(vapos+runningins-runningdel)]
						readprevious=seq[(vapos+runningins-runningdel-8):(vapos+runningins-runningdel-1)]
						if genomeprevious==readprevious:
							#print('-->2variant: reference %s alternate %s, read %s,  genomeprevious %s, readprevious %s, read-coord-of-var is %s, pointer %s, running deletions is %s and running insertions is %s'%(varef,vaalt, readvar, genomeprevious, readprevious, (vapos+runningins-runningdel) ,pointer, runningdel, runningins))
							found+=1
							if readvar==varef:
								#print('reference')
								call_list.append('ref')
								refcounter+=1
							elif readvar==vaalt:
								#print('alt')
								call_list.append('alt')
								altcounter+=1
							else:
								errorcounter+=1
								#call_list.append('error')
						else:
							oldvar=readvar
							plusone=seq[(vapos+runningins-runningdel-7):(vapos+runningins-runningdel)]
							plustwo=seq[(vapos+runningins-runningdel-6):(vapos+runningins-runningdel+1)]
							plusthree=seq[(vapos+runningins-runningdel-5):(vapos+runningins-runningdel+2)]
							minusone=seq[(vapos+runningins-runningdel-9):(vapos+runningins-runningdel-2)]
							minustwo=seq[(vapos+runningins-runningdel-10):(vapos+runningins-runningdel-3)]
							minusthree=seq[(vapos+runningins-runningdel-11):(vapos+runningins-runningdel-4)]
							if genomeprevious==plusone:
								context=seq[(vapos+runningins-runningdel-7):(vapos+runningins-runningdel+1)]
							elif genomeprevious==plustwo:
								context=seq[(vapos+runningins-runningdel-6):(vapos+runningins-runningdel+2)]
							elif genomeprevious==plusthree:
								context=seq[(vapos+runningins-runningdel-5):(vapos+runningins-runningdel+3)]
							elif genomeprevious==minusone:
								context=seq[(vapos+runningins-runningdel-9):(vapos+runningins-runningdel-1)]
							elif genomeprevious==minustwo:
								context=seq[(vapos+runningins-runningdel-10):(vapos+runningins-runningdel-2)]
							elif genomeprevious==minusthree:
								context=seq[(vapos+runningins-runningdel-11):(vapos+runningins-runningdel-3)]
							else:
								match=seq.find(longprevious)
								if match != -1:
									context=seq[match:match+13]
								else:
									readvar='null'
							if readvar !='null':
								readvar=context[::-1][0]#get last item
								#print('-->2variant: reference %s alternate %s, read %s,  genomeprevious %s, readprevious %s, read-coord-of-var is %s, pointer %s, running deletions is %s and running insertions is %s'%(varef,vaalt, readvar, genomeprevious, readprevious, (vapos+runningins-runningdel) ,pointer, runningdel, runningins))
								found+=1
								if readvar==varef:
									#print('reference')
									call_list.append('ref')
									refcounter+=1
								elif readvar==vaalt:
									#print('alt')
									call_list.append('alt')
									altcounter+=1
								else:
									errorcounter+=1
									#call_list.append('error')

				#print('readname %s,good align, %s %s, and length of %s,  %s variants in this region of which found %s, with %s reference and %s alternate and %s errors (match neither alt nor ref)'%(readname,mapcontig,mappos,len(seq),vc,found,refcounter,altcounter, errorcounter))
				if (refcounter+altcounter)!=0:
					fracref=float((refcounter/(refcounter+altcounter)))
					refFlag=0
					altFlag=0
					recombFlag=0
					if fracref > 0.75:
						refFlag=1
					else:
						altFlag=1
				
					if len(call_list)>20:
						first10=[i for i in call_list[0:10]]
						last10=[i for i in call_list[len(call_list)-10:len(call_list)]]
						firstset=set(first10)
						lastset=set(last10)
						if len(firstset)==1 and len(lastset)==1:#ensure no switching in the first or last 10
							if firstset!=lastset:#first and last do not match
								recombFlag=1
								recombcounter+=1
					callHandle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(readname,mapcontig,mappos,refFlag,altFlag,recombFlag,refcounter,altcounter,errorcounter,(len(seq)+mappos),call_list))

		if x%1000==0:
			print('%s thousand reads processed and %s recombinant so far. Currently on contig %s'%(x/1000,recombcounter,mapcontig))
print('of %s total reads %s were good and long enough'%(x,m))
samHandle.close()
callHandle.close()
