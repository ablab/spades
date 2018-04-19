import sys
import collections
import numpy as np
from sklearn.naive_bayes import MultinomialNB

# creates scikit_learn classifier based on input plasmid and chromosomal tblouts 


plasmid_file =   "/Nancy/mrayko/plasmid_prediction_test/HMM_refseq/all_refseq_plasmids_PfamA_001_tblout"   #sys.argv[1] "plasmids_1000.tblout" 
chromosome_file = "/Nancy/mrayko/plasmid_prediction_test/downsampling_20perc/chrom_10percen/total_filtered.fna_tblout"  #sys.argv[2]  "chromosomes_1000.tblout"

with open("pfam_names.list", "r") as infile:
        pfam_list=infile.readlines()
pfam_list = [i.strip() for i in pfam_list] 




def get_annot_from_tblout(tblout_pfam):
# Takes tblout file, returns list of contig names + hmms (e.g. ['NODE_4_length_6662_cov_295.444039_cutoff_5', 'ABC2_membrane', 'HTH_3', 'Rep_3', 'Relaxase', 'MobC'])    
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()

    tblout_pfam = [i.split() for i in tblout_pfam] 
    contigs = []
    for i in tblout_pfam[3:-10]:
        name = [i[2].rsplit("_", maxsplit=1)[0]]
        if name not in contigs: 
        	contigs.append(name)
    
    
    for item in contigs: # In table - add new field, and for each row
        pr_list=[]
        for j in tblout_pfam[3:-10]:  #check tblout (last 10 - service strings)
            if j[2] and item[0]: # if there is a protein in tblout and SRR in table
                if j[2] not in pr_list: #  if we didn't see this protein yet
                    if item[0] in j[2]  and float(j[4]) < 0.01  : # if there's the SRR in protein ID and e-value < 0.01
                        pr_list.append(j[2]) # add this protein to the list for given SRR
                        item.append(str(j[0])) # add to this new field
    
    return (contigs)


def get_annot_from_tblout(tblout_pfam):
# Takes tblout file, returns list of contig names + hmms (e.g. ['NODE_4_length_6662_cov_295.444039_cutoff_5', 'ABC2_membrane', 'HTH_3', 'Rep_3', 'Relaxase', 'MobC'])    
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()

    tblout_pfam = [i.split() for i in tblout_pfam] 
    contigs = collections.OrderedDict()
    for i in range(3, len(tblout_pfam)-10):
      #  print (tblout_pfam[i])
        name = tblout_pfam[i][2].rsplit("_", maxsplit=1)[0]  # get name
        if name not in contigs: 
            contigs[name]=[tblout_pfam[i][0]]
        else:
            if tblout_pfam[i][2] != tblout_pfam[i-1][2]:
                contigs[name].append(tblout_pfam[i][0])

    return contigs
            





def create_vector_pfams(hmms): # list of hmm lists
   vector=[len(pfam_list)*[0]]*len(hmms)
   for i in range(0,len(hmms)): # take each hit
       for j in hmms[i]:
       	#print (j)
        hit_index = pfam_list.index(j)
        vector[i][hit_index]+=1
   return vector


    
print ("Extracting hmms from tblout file...")

plasmid_train = get_annot_from_tblout(plasmid_file)
chromosome_train = get_annot_from_tblout(chromosome_file)

print (len(plasmid_train))
print (len(chromosome_train))


# ok, let's assume we have two training sets. How to train classifier?
#Combine both datasets:

train=collections.OrderedDict(plasmid_train, **chromosome_train)

#print(train)

print ("Initializing matrix...")


X=[len(pfam_list)*[0]]*len(train)  # initialize matirx of n rows of zeroes, len = num of pfams

print ("Populating matrix...")

counter=0
for key, value in train.items(): # take each sample  key, value in d.items():
#	print (key)
#	print (value)
	for j in value: # take each hit
	    hit_index = pfam_list.index(j)  # int(str([x for x in range(len(pfam_list)) if pfam_list[x]==j]))   #int(i for i,x in enumerate(pfam_list) if x == j) # pfam_list.index(j)  
	    X[counter][hit_index]+=1
	counter+=1

y=["Plasmid"]*len(plasmid_train)+["Chromosome"]*len(chromosome_train)



print ("MultinomialNB training...")

# Train classificator
clf = MultinomialNB()
clf.fit(X, y)

# Test
a=create_vector_pfams([['ABC2_membrane', 'HTH_3', 'Rep_3', 'Relaxase', 'MobC']])

print(clf.predict(a))
print(clf.predict_proba(a))




import pickle
# save the classifier
with open('my_dumped_classifier.pkl', 'wb') as fid:
    pickle.dump(clf, fid)    
