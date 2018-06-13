import sys
from operator import itemgetter

with open("/Nancy/mrayko/db/pfam/Pfam-A.hmm", 'r') as infile:
    hmms=infile.readlines()
    hmms = [i.strip() for i in hmms] 

hmm_nc={}
for i in range (len(hmms)):
    if hmms[i][:4] == "NAME":
        hmm_nc[hmms[i].split()[1]]=hmms[i+16].split()[1]


t=sys.argv[1]
#t="final_contigs_cat.fasta_tblout"
t1="/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/best_k/all_hmms_sorted"


with open(t) as f:
    tblout_pfam = f.read().splitlines()

with open(t1) as f1:
    hmms = f1.read().splitlines()


hmms = [i.split("\t") for i in hmms]

pl_list=[]
for i in hmms:
	pl_list.append(i[0])

tblout_pfam = [i.split() for i in tblout_pfam] 

def get_plasmids_number(pl_hmm_list):
  plasmids=set()
  pls=set(pl_hmm_list)
  for i in tblout_pfam[3:-10]: #watch out for tblout structure
    if float (i[5]) >= float(hmm_nc[i[2]]):
    	if i[2] in pls:
    	    plasmids.add(i[0].rsplit('_', 1)[0])
  return len(plasmids)


#print (get_plasmids_number(pl_list[:377]))
#print (get_plasmids_number(pl_list[:378]))
#print (get_plasmids_number(pl_list[:379]))

k=[154, 276, 378, 496, 570, 698]


for i in k: #range (0, 1641, 10):   #len(pl_list),500):

	print ("hmms: " + str(i) +  " Plasmids: " + str(get_plasmids_number(pl_list[:i])) +  " K: "+str(hmms[i-1][5]))



