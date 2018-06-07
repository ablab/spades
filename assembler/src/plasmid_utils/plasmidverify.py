import os, errno
import sys
import ntpath
import argparse
import subprocess
from joblib import Parallel, delayed
from glob import glob
from Bio import SeqIO
from parse_blast_xml import parser
from classifier import naive_bayes
from classifier import scikit_multNB

base = os.path.basename(sys.argv[1])
name_file = os.path.splitext(base)[0]

outdir = sys.argv[2]

try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

name=os.path.join(outdir, name_file)

#hmm= "/Nancy/mrayko/db/pfam/Pfam-A.hmm"
#hmm="/Nancy/mrayko/db/plasmid_specific_pfam/1641_hmm_list.hmm" # up to k=1
hmm= "/Nancy/mrayko/db/plasmid_specific_pfam/378_10fold_plasmid_HMMs.hmm"
list378="/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/378_hmms.txt" 
blastdb=" /Bmo/ncbi_nt_database/nt"


hmmsearch=" /Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch"
prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"
cbar="/Nancy/mrayko/Libs/cBar.1.2/cBar.pl"


# run hmm
os.system (prodigal + " -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
os.system (hmmsearch + " --noali --cut_nc  -o "+name+"_out_pfam --tblout "+name+"_tblout --cpu 10 "+ hmm + " "+name+"_proteins.fa")
os.system ("tail -n +4 " + name +"_tblout | head -n -10 | awk '{print $1}'| sed 's/_[^_]*$//g'| sort | uniq > " + name +"_plasmid_contigs_names.txt")


# parse hmms
with open(list378, "r") as infile1:
        list378=infile1.readlines()

with open(name + "_tblout", "r") as infile2:
        tblout_pfam=infile2.readlines()

list378 = [i.strip() for i in list378] 
for i in tblout_pfam:
    i=[x for x in i if x]


# Create table backbone.
table =[]

ids=[]
records = list(SeqIO.parse(sys.argv[1], "fasta"))
for i in records: 
    ids.append(i.id)  # take all fasta ids and append to table

for i in ids:
  table.append(i.split())


tblout_pfam = [i.split() for i in tblout_pfam] 
# add list of domains
for item in table: # In table - add new field, and for each row
    item.append([])
    for j in tblout_pfam[3:-10]:  #check tblout (last 10 - service strings)
        if j[0] and item[0]: # if there is a protein in tblout and SRR in table
#            if j[2] not in pr_list: #  if we didn't see this protein yet
                if item[0] == j[0].rsplit('_', 1)[0]: # if there's the SRR in protein ID and bitscore > 20
                    item[-1].append(str(j[2])) # add to this new field
    # add our list378
    item.insert(-1,[])
    for i in item[-1]:
       if i in list378: 
           item[-2].append(i)
    if len(item[-2])==0: 
        item[-2]=("-")
    if len(item[-1])==0:
        item[-1]=("-")

    item[-1]=' '.join(item[-1]) # turn list of genes into the string
    item[-2]=' '.join(item[-2])


# Add HMM results
for i in table:
#    print(i[1])
    if i[1]!="-":
        i.append("hmm+")
    else:
        i.append("hmm-")

# Output

with open(name + "_result_table.tsv", "w") as outfile:
    for i in table:
         outfile.write('\t'.join(i)+"\n")
