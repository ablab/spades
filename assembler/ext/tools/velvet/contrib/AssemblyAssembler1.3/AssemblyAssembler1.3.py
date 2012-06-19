#!/usr/bin/python
#
#       AssemblyAssembler.pl
#
#		Conduct velvet assemblies across a range of kmer values, find parameter region 
#		with good assembly stats, conduct additional assemblies in that region and then
#		assemble contigs from all previous assemblies in one final assembly.  Potentially 
#		useful for de novo transcriptome assembly.  
#
#       Copyright 2010 Jacob Crawford <jc598@cornell.edu>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#
# ==========================================================================================

## import modules
import sys
import os

## Checking version of python
if sys.version_info < (2, 5):
    sys.exit('\n Incorrect version of Python.  Please update to 2.5 or greater \n')

## Setting working directory
parDir = os.getcwd()

print str('============================================== \n')
## Retrieve user input
ptester = [0,0,0,0,0,0,0,0,0]
increads = str('n')
print '\nRetrieving your input values \n'
for i in  range(len(sys.argv)):
	if (sys.argv[i] == '-s'):
		if(int(sys.argv[i+1])%2 == 1):
			hashs = sys.argv[i+1]
			ptester[0] = 1
		else:
			sys.exit('Invalid hash value.  Please enter odd number')
	elif (sys.argv[i] == '-e'):
		if(int(sys.argv[i+1])%2 == 1):
			hashe = sys.argv[i+1]
			ptester[1] = 1
		else:
			sys.exit('Invalid hash value.  Please enter odd number')
	elif (sys.argv[i] == '-f'):
		## Velvethcall must be entered as a string at the command line with ' ' or " "
		velvethcall = sys.argv[i+1]
		ptester[2] = 1
	elif (sys.argv[i] == '-v'):
		velvDir = sys.argv[i+1]
		ptester[3] = 1
	elif (sys.argv[i] == '-c'):
		covcut = sys.argv[i+1]
		ptester[4] = 1
	elif (sys.argv[i] == '-t'):
		thresh = sys.argv[i+1]
		ptester[5] = 1
	elif (sys.argv[i] == '-i'):
		inslgth = sys.argv[i+1]
		ptester[6] = 1
	elif (sys.argv[i] == '-m'):
		expcov = sys.argv[i+1]
		ptester[7] = 1
	elif (sys.argv[i] == '-r'):
		increads = sys.argv[i+1].lower()
		ptester[8] = 1

## Check and set hash value range
print 'Checking your hash value range...'
if ((int(hashe)-int(hashs))>=16):
	print 'Looks good\n'
	hashrange = range(int(hashs),(int(hashe)+1),2)
	hashlist = [hashrange[0],hashrange[int(round(len(hashrange)*0.12))],hashrange[int(round(len(hashrange)*0.25))],hashrange[int(round(len(hashrange)*0.37))],hashrange[int(round(len(hashrange)*0.5))],hashrange[int(round(len(hashrange)*0.6))],hashrange[int(round(len(hashrange)*0.7))],hashrange[int(round(len(hashrange)*0.82))],hashrange[len(hashrange)-1]]
else:
	sys.exit('Insufficient hash range. This routine is only useful if kmer range >= 16. \n\n')
		
## Checking to make sure minimum parameters are in place
print 'Checking to make sure minimum parameters are in place...'
if (sum(ptester[0:4]) < 4):
	sys.exit('Incorrect command line entry.  Please enter parameters values for -s, -e, -f -v and -d \n')
elif(sum(ptester[0:4]) == 4): 
	if(sum(ptester[4:6])==0):
		print str('Okay, will conduct assemblies across range of kmer values from '+hashs+' to '+hashe+' with no coverage cutoff \n')
	elif(sum(ptester[4:6])==1):
		sys.exit('Incorrect command line entry.  Not quite enough information. Please enter both a coverage cutoff value (or auto) and the lowest kmer value to implement the cutoff \n')
	elif(sum(ptester[4:6])==2):
		print str('Okay, will conduct assemblies across range of kmer values from '+hashs+' to '+hashe+' with coverage cutoff of '+covcut+' implemented at kmer values between '+thresh+' and '+hashe+ '\n') 
	if(sum(ptester[6:8])==0):
		print str('Running assemblies with unpaired reads \n')	
	elif(sum(ptester[6:8])==1):
		if(ptester[7]==0):
			sys.exit('Incorrect command line entry.  Not quite enough information. Please enter both a expected coverage value (or auto) and the expected length of the total sequenced fragment (see Velvet manual) \n')	
	elif(sum(ptester[6:8])==2):
		print str('Will run assemblies with paired reads \n')	


print str('==============================================')
print str('Step 1: Exploratory Assemblies \n')
## Constructing Velvetg call
vgcallb = str(velvDir+'/velvetg ./')
if(sum(ptester[6:8])==2):
	## use paired reads
	vgcall = str(vgcallb+' -ins_length '+inslgth+' -exp_cov '+expcov)
elif(sum(ptester[6:8])<2):
	## use un-paired reads
	vgcall = vgcallb
if(sum(ptester[4:6])==2):
	## Call below threshold
	vgcall = vgcall
	## Call above threshold
	vgcallat = str(vgcall+' -cov_cutoff '+covcut)

## Final Velvetg calls	
vgcallat = str(vgcall+' -unused_reads yes >> '+parDir+'/FinalDir/GrandVelvetLog.txt')
vgcallatOpt = str(vgcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt')
vgcall = str(vgcall+' -unused_reads yes >> '+parDir+'/FinalDir/GrandVelvetLog.txt')
vgcallOpt = str(vgcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt')

	
## Make directories, one for each hash value
try:
	os.system(str('rm -rf SceneOfTheCrime'))
	os.system(str('rm -rf FinalDir'))
	os.mkdir('SceneOfTheCrime')
	os.mkdir('FinalDir')
	os.mkdir('FinalDir/UU')
	os.chdir(str(parDir+'/SceneOfTheCrime'))
except:
	os.mkdir('SceneOfTheCrime')
	os.mkdir('FinalDir')
	os.mkdir('FinalDir/UU')
	os.chdir(str(parDir+'/SceneOfTheCrime'))	

hashes = [hashlist[0],hashlist[2],hashlist[4],hashlist[6],hashlist[8]]

## Call velveth and velvetg for hash values across coarse sampling of range and determining best assembly
print 'Running velvet with intitial set of hash values \n'
contigs = []
unu = []
MAXCONT = {}
maxcheck = [0,0,0,0,0]
for i in range(5):
	os.mkdir(str('Dir'+str(hashes[i])))
	os.chdir(str(os.getcwd()+'/Dir'+str(hashes[i])))
	print 'Running velveth with k= '+str(hashes[i])
	os.system(str(velvDir+'/velveth ./ '+str(hashes[i])+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	if(sum(ptester[4:6])==0):
		print 'Running velvetg on k= '+str(hashes[i])+' velveth files'
		os.system(vgcall)
	elif(sum(ptester[4:6])==2):
		if(hashes[i] < thresh):
			print 'Running velvetg on k= '+str(hashes[i])+' velveth files'
			os.system(vgcall)
		elif(hashes[i] >= thresh):
			print 'Running velvetg with coverage cutoff on k= '+str(hashes[i])+' velveth files'
			os.system(vgcallat)
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(hashes[i])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(hashes[i])))
	contigs.append(str('contigs'+str(hashes[i])+'.fa '))
	if(int(i)%2 == 0):
		os.system(str('mv UnusedReads.fa '+parDir+'/FinalDir/UU/unused'+str(hashes[i])+'.fa'))
		unu.append(str('./unused'+str(hashes[i])+'.fa '))	
	elif(int(i) == 1):
		os.system('rm UnusedReads.fa')
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(hashes[i])))
	input = open(str(parDir+'/FinalDir/Log'+str(hashes[i])),'r')
	for j in input:
		j = j.split(' ')
		if (j[0] == 'Final'):
			if not MAXCONT.has_key(hashes[i]):
				MAXCONT[hashes[i]] = 0
			MAXCONT[hashes[i]] = int(j[10][:-1])
			maxcheck[i] = 1
	input.close()

if (sum(maxcheck) < 5):
	sys.exit('One or more of the exploratory assemblies failed.  Please examine log files to determine the cause and try again.')

print str('==============================================')
print str('Step 2: Optimisation Assemblies \n')
## Make directories and do additional velvet runs
BEST = max(MAXCONT, key = MAXCONT.get)
print str('Of the first 5 assemblies k='+str(BEST)+' had the longest contig')
print str('Sampling deeper around '+str(BEST)+'\n')
if (BEST == hashlist[0]):
	print 'Running velveth with k= '+str(hashlist[1])
	hashes.append(hashlist[1])
	os.mkdir(str('Dir'+str(hashlist[1])))
	os.chdir(str(os.getcwd()+'/Dir'+str(hashlist[1])))
	os.system(str(velvDir+'/velveth ./ '+str(hashlist[1])+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(hashlist[1])+' velveth files'
	if(sum(ptester[4:6])==0):
		os.system(vgcallOpt)
	elif(sum(ptester[4:6])==2):
		if(BEST < thresh):
			os.system(vgcallOpt)
		elif(BEST >= thresh):
			os.system(vgcallatOpt)
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(hashlist[1])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(hashlist[1])))
	contigs.append(str('contigs'+str(hashlist[1])+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(hashlist[1])))
	
elif (BEST == hashlist[8]):
	print 'Running velveth with k= '+str(hashlist[7])
	hashes.append(hashlist[7])
	os.mkdir(str('Dir'+str(hashlist[7])))
	os.chdir(str(os.getcwd()+'/Dir'+str(hashlist[7])))
	os.system(str(velvDir+'/velveth ./ '+str(hashlist[7])+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(hashlist[7])+' velveth files'
	if(sum(ptester[4:6])==2):
		if(BEST < thresh):
			os.system(vgcallOpt)
		elif(BEST >= thresh):
			os.system(vgcallatOpt)
	elif(sum(ptester[4:6])==0):
		os.system(vgcallOpt)
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(hashlist[7])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(hashlist[7])))
	contigs.append(str('contigs'+str(hashlist[7])+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(hashlist[7])))
	
else:
	for i in range(9):
		if (hashlist[i] == BEST):
			HASHUP = hashlist[i+1]
			HASHDWN = hashlist[i-1]		
	hashes.append(HASHUP)
	hashes.append(HASHDWN)
	os.mkdir(str('Dir'+str(HASHUP)))
	os.chdir(str(os.getcwd()+'/Dir'+str(HASHUP)))
	print 'Running velveth with k= '+str(HASHUP)
	os.system(str(velvDir+'/velveth ./ '+str(HASHUP)+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(HASHUP)+' velveth files'
	if(sum(ptester[4:6])==0):
		os.system(vgcallOpt)
	elif(sum(ptester[4:6])==2):
		if(HASHUP < thresh):
			os.system(vgcallOpt)
		elif(HASHUP >= thresh):
			os.system(vgcallatOpt)
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(HASHUP)+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(HASHUP)))
	contigs.append(str('./contigs'+str(HASHUP)+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(HASHUP)))
	os.mkdir(str('Dir'+str(HASHDWN)))
	os.chdir(str(os.getcwd()+'/Dir'+str(HASHDWN)))
	print 'Running velveth with k= '+str(HASHDWN)
	os.system(str(velvDir+'/velveth ./ '+str(HASHDWN)+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(HASHDWN)+' velveth files'
	if(sum(ptester[4:6])==0):
		os.system(vgcallOpt)
	elif(sum(ptester[4:6])==2):
		if(HASHDWN < thresh):
			os.system(vgcallOpt)
		elif(HASHDWN >= thresh):
			os.system(vgcallatOpt)
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigs'+str(HASHDWN)+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/Log'+str(HASHDWN)))
	contigs.append(str('contigs'+str(HASHDWN)+'.fa '))
	os.chdir('..')
	os.system(str('rm -rf Dir'+str(HASHDWN)))

## Assemble unused reads 	
os.chdir(str(parDir+'/FinalDir/UU'))
print '\nConducting assembly of Unused reads from exploratory assemblies'
print 'Running velveth with k= '+str(hashlist[4])
os.system(str(velvDir+'/velveth ./ '+str(hashlist[4])+'-fasta -short '+''.join(unu)+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
print 'Running velvetg on k= '+str(hashlist[4])+' velveth files'
os.system(str(velvDir+'/velvetg ./ >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
os.system(str('mv contigs.fa '+parDir+'/FinalDir/UnusedOut.fa'))
os.system(str('mv Log '+parDir+'/FinalDir/LogUU'))
os.chdir('..')
os.system('rm -rf UU')

print str('==============================================')
print str('Step 3: Summary Assemblies \n')
## Do first round of summary assemblies
sumhashes = [hashlist[2],hashlist[4],hashlist[6]]
contigsSum = []
for i in range(3):
	os.mkdir(str(parDir+'/FinalDir/Sum'))
	os.chdir(str(parDir+'/FinalDir/Sum'))
	files = str('')
	for j in range(len(contigs)):
		files = str(files+parDir+'/FinalDir/'+contigs[j]+' ')
	files = str(files+parDir+'/FinalDir/UnusedOut.fa')
	print 'Running velveth with k= '+str(sumhashes[i])
	if(increads==str('n')):
		os.system(str(velvDir+'/velveth ./ '+str(sumhashes[i])+' -fasta -long '+files+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
		print 'Running velvetg on k= '+str(sumhashes[i])+' velveth files'
		os.system(str(velvDir+'/velvetg ./ >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	elif(increads==str('y')):
		os.system(str(velvDir+'/velveth ./ '+str(sumhashes[i])+' -fasta -long '+files+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
		print 'Running velvetg on k= '+str(sumhashes[i])+' velveth files'
		if(sum(ptester[6:8])<2):
			os.system(str(velvDir+'/velvetg ./ >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))	
		elif(sum(ptester[6:8])==2):
			os.system(str(velvDir+'/velvetg ./ -ins_length '+inslgth+' -exp_cov '+expcov+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	os.system(str('mv contigs.fa '+parDir+'/FinalDir/contigsSum'+str(sumhashes[i])+'.fa'))
	os.system(str('mv Log '+parDir+'/FinalDir/LogSum'+str(sumhashes[i])))
	contigsSum.append(str('./contigsSum'+str(sumhashes[i])+'.fa'))
	os.chdir('..')
	os.system('rm -rf Sum')

## Do final summary assembly
os.chdir(str(parDir+'/FinalDir'))
print '\nConducting final summary assembly'
print 'Running velveth with k= '+str(hashlist[4])
if(increads==str('n')):
	os.system(str(velvDir+'/velveth ./ '+str(hashlist[4])+' -fasta -long '+' '.join(contigsSum)+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(hashlist[4])+' velveth files'
	os.system(str(velvDir+'/velvetg ./ >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
elif(increads==str('y')):
	os.system(str(velvDir+'/velveth ./ '+str(hashlist[4])+' -fasta -long '+' '.join(contigsSum)+' '+velvethcall+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
	print 'Running velvetg on k= '+str(hashlist[4])+' velveth files'
	if(sum(ptester[6:8])<2):
		os.system(str(velvDir+'/velvetg ./ >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))	
	elif(sum(ptester[6:8])==2):
		os.system(str(velvDir+'/velvetg ./ -ins_length '+inslgth+' -exp_cov '+expcov+' >> '+parDir+'/FinalDir/GrandVelvetLog.txt'))
os.system('mv contigs.fa Finalcontigs.fa')
os.system('mv Log LogSumFinal')
	
## Remove interim directories and their contents
os.chdir('..')
os.system('rm -rf SceneOfTheCrime')

print str('============================================== \n')	
print 'Assembly Assembler Job Complete \n' 

input = open(str(parDir+'/FinalDir/LogSumFinal'),'r')
for i in input:
	i = i.split(' ')
	if (i[0] == 'Final'):
		print ' '.join(i)
input.close()
























