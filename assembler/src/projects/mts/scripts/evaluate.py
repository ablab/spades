#!/usr/bin/env python2
#Evaluation script from MyCC

import os
import os.path
import re
import numpy as np
import sys
import commands
argv=sys.argv
split=False
coords=''
if '-h' in argv or len(argv)==1:
    print 'Please input two files containing reference assignment and assembly, followed by target name (e.g., fasta).'
    print 'Usage:'
    print 'Evaluate.py [reference assignment] [assembly] [target name] [options]'
    print 'option:'
    print '-split To split header by space.'
    print '-plot To plot references based on the MyCluster.coords.'
    print '-h Help.'
    exit()
if '-split' in argv:
    split=True
if '-plot' in argv:
    coords='MyCluster.coords'
myFile=argv[1]
ctgfile=argv[2]
target=argv[3]
#myFile='./Species.txt'
#target='bin'
outP='Precision'
outS='Sensitivity'
ClusterSeq=dict()
Seqs=[]
SeqLen=dict()
myfilelist=[]
for i in os.listdir('./'):
    if target in i:
        bases=0
        sample,_=os.path.splitext(i) #Hack for MTS compatibility
        f=open(i)
        for j in f:
            if '>' in j:
                h=sample+"-"+j.strip().replace('>','')
                if split:
                    h=h.split()[0]
                SeqLen[h]=0
                continue
            bases+=len(j.strip())
            SeqLen[h]+=len(j.strip())
        f.close()
        #if bases < 200000:
        #    print 'skip: %s %d'%(i,bases)
        #    continue
        myfilelist.append(i)
if not myfilelist:
    print 'Please make sure your current directory is YYYYMMDD_MMSS_mer_lt for Evaluate.py.'
    exit()

for i in myfilelist:
    ClusterSeq[i]=[]
    sample,_=os.path.splitext(i) #Hack for MTS compatibility
    f=open(i)
    for j in f:
        if '>' in j:
            j=j.strip().replace('>','')
            if split:
                j=j.split()[0]
            ClusterSeq[i].append(sample+"-"+j)
            Seqs.append(j)
    f.close()

CtgInfo2Ref=dict()
f=open(myFile)
#f.readline()
for i in f:
    [ctginfo,ref]=i.strip().split('\t')
    if split:
        ctginfo=ctginfo.split()[0]
    CtgInfo2Ref[ctginfo]=ref
f.close()

AllCtgSeqLen=dict()
f=open(ctgfile)
for i in f:
    i=i.strip()
    if '>' in i:
        h=i.replace('>','')
        if split:
            h=h.split()[0]
        AllCtgSeqLen[h]=0
        continue
    AllCtgSeqLen[h]+=len(i)
f.close()

#Precision
Precision=dict()
Sentivity=dict()
Precisionw=dict()
Sentivityw=dict()
totalSeqs=0
for cluster,cSeqs in ClusterSeq.items():
    for s in cSeqs:
        #ss=re.findall('(\[.+\])',s)[0].replace('[','').replace(']','')
        try:
            ref=CtgInfo2Ref[s]
            totalSeqs+=1
        except:
            #print cluster,s
            continue
        myLen=SeqLen[s]
        Precision.setdefault(cluster,[]).append(ref)
        Sentivity.setdefault(ref,[]).append(cluster)
        Precisionw.setdefault(cluster,[]).extend([ref]*myLen)
        Sentivityw.setdefault(ref,[]).extend([cluster]*myLen)

AllCtgSeqLenWithRef=dict()
for h,v in AllCtgSeqLen.items():
    if CtgInfo2Ref.get(h):
        AllCtgSeqLenWithRef[h]=v

AllCtgSeq=len(AllCtgSeqLenWithRef.keys())
AllCtgSeqWeight=sum(AllCtgSeqLenWithRef.values())

print 'No. of reference genomes: %d'%len(Sentivity.keys())
print 'No. of bins in evaluation: %d'%len(Precision.keys())
print 'No. of sequences assigned reference: %d'%AllCtgSeq
print 'No. of binned sequences: %d'%len(SeqLen.keys())


fw=open(outP+'.txt','w')
Precision1_up=[]
Precision1_down=[]
for cluster,refs in Precision.items():
    toprefCount=sorted([refs.count(x) for x in set(refs)],reverse=True)[0]
    myp=toprefCount/float(len(refs))
    Precision1_up.append(toprefCount)
    Precision1_down.append(float(len(refs)))
    fw.write(cluster+':%f(%d/%d)'%(myp,toprefCount,len(refs))+'\t'+('\t'.join(['%s:%d'%(x,refs.count(x)) for x in set(refs)]))+'\n')
fw.close()
fw=open(outP+'.w.txt','w')
Precisionw1_up=[]
Precisionw1_down=[]
for cluster,refs in Precisionw.items():
    toprefCount=sorted([refs.count(x) for x in set(refs)],reverse=True)[0]
    myp=toprefCount/float(len(refs))
    Precisionw1_up.append(toprefCount)
    Precisionw1_down.append(float(len(refs)))
    fw.write(cluster+':%f(%d/%d)'%(myp,toprefCount,len(refs))+'\t'+('\t'.join(['%s:%d'%(x,refs.count(x)) for x in set(refs)]))+'\n')
fw.close()
print 'Precision: %f, %f'%(sum(Precision1_up)/float(sum(Precision1_down)),sum(Precisionw1_up)/float(sum(Precisionw1_down)))

fw=open(outS+'.txt','w')
Sentivity1_up=[]
Sentivity1_down=[]
for ref,clusters in Sentivity.items():
    topClusterCount=sorted([clusters.count(x) for x in set(clusters)],reverse=True)[0]
    mys=topClusterCount/float(len(clusters))
    Sentivity1_up.append(topClusterCount)
    Sentivity1_down.append(float(len(clusters)))
    fw.write(ref+':%f(%d/%d)'%(mys,topClusterCount,len(clusters))+'\t'+('\t'.join(['%s:%d'%(x,clusters.count(x)) for x in set(clusters)]))+'\n')
fw.close()
fw=open(outS+'.w.txt','w')
Sentivityw1_up=[]
Sentivityw1_down=[]
for ref,clusters in Sentivityw.items():
    topClusterCount=sorted([clusters.count(x) for x in set(clusters)],reverse=True)[0]
    mys=topClusterCount/float(len(clusters))
    Sentivityw1_up.append(topClusterCount)
    Sentivityw1_down.append(float(len(clusters)))
    fw.write(ref+':%f(%d/%d)'%(mys,topClusterCount,len(clusters))+'\t'+('\t'.join(['%s:%d'%(x,clusters.count(x)) for x in set(clusters)]))+'\n')
fw.close()
#print 'Sensitivity: %f, %f'%(sum(Sentivity1_up)/float(sum(Sentivity1_down)),sum(Sentivityw1_up)/float(sum(Sentivityw1_down)))
print 'Sensitivity: %f, %f'%(sum(Sentivity1_up)/float(AllCtgSeq),sum(Sentivityw1_up)/float(AllCtgSeqWeight))


if coords:
    try:
        f=open(coords)
    except:
        print 'Please make sure your current directory is YYYYMMDD_MMSS_mer_lt for Evaluate.py.'
        exit()
    fw=open('MyReference.coords','w')
    for i in f:
        #[h,x,y,c]=i.strip().split()
        tmp=i.strip().split('\t')
        h=tmp[0];x=tmp[1];y=tmp[2];c=tmp[3]
        if split:
            h=h.split()[0]
        try:
            ref=CtgInfo2Ref[h]
            fw.write('%s\t%s\t%s\t%s\n'%(h,x,y,ref))
        except:
            continue
    fw.close()
    f.close()

    comm='getCentroid.py MyReference.coords'
    (s,o)=commands.getstatusoutput(comm)
    comm='ToDrawPDF.py MyReference.coords MyReference.pdf'
    (s,o)=commands.getstatusoutput(comm)
