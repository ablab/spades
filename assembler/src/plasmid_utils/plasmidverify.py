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

hmm= "hmms/378_10fold_plasmid_HMMs.hmm"
list378="hmms/378_hmms.txt" 


#hmmsearch="/Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch"
#prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"
#cbar="/Nancy/mrayko/Libs/cBar.1.2/cBar.pl"


# run hmm
os.system ("prodigal -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
os.system ("hmmsearch --noali --cut_nc  -o "+name+"_out_pfam --tblout "+name+"_tblout --cpu 10 "+ hmm + " "+name+"_proteins.fa")
os.system ("tail -n +4 " + name +"_tblout | head -n -10 | awk '{print $1}'| sed 's/_[^_]*$//g'| sort | uniq > " + name +"_plasmid_contigs_names.txt")


# parse hmms
with open(list378, "r") as infile1:
        list378=infile1.readlines()

with open(name + "_tblout", "r") as infile2:
        tblout_pfam=infile2.readlines()

list378 = [i.strip() for i in list378] 
set378 = set(list378)

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



# read file one time, create dict - plasmid and all genes that fit.
tblout_pfam = [i.split() for i in tblout_pfam] 

plasmid_hits={}

for i in tblout_pfam[3:-10]:
    if i[0] and i[2]:
        plasmid_hits.setdefault(i[0].rsplit('_', 1)[0],[]).append(i[2]) 


# Add genes for each contig
for item in table: 
    genes=list(plasmid_hits.get(item[0], ""))
    hits=[]
    for i in genes:
        if i in set378:
            hits.append(i)

    if len(hits)==0: 
        hits="-"
    if len(genes)==0:
        genes="-"

    item.append(' '.join(hits))
    item.append(' '.join(genes))
    if hits!="-":
        item.append ("hmm+")
    else: 
        item.append ("hmm-")

# Output

with open(name + "_result_table.tsv", "w") as outfile:
    for i in table:
         outfile.write('\t'.join(i)+"\n")
