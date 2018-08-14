#!/usr/bin/env python

import os, errno
import sys
import argparse

from Bio import SeqIO
from parse_blast_xml import parser


import ntpath
import subprocess
from joblib import Parallel, delayed
from glob import glob
from Bio import SeqIO
from classifier import naive_bayes
from classifier import scikit_multNB



def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="HMM-based plasmid verification script")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('-f', required = True, help='Input fasta file')
    parser.add_argument('-o', required = True, help='Output directory')
    parser.add_argument('-b', help='Run BLAST on input contigs', action='store_true')
    parser.add_argument('-c', help='Run Scikit-Learn Classifier on input contigs', action='store_true')
    parser.add_argument('-db', help='Path to BLAST db')
    parser.add_argument('-hmm', help='Path to plasmid-specific HMMs')    
    parser.add_argument('-K', help='K for plasmid/chromosome ratio, K>10') 
    return parser.parse_args()



args = parse_args(sys.argv[1:])


base = os.path.basename(args.f)
name_file = os.path.splitext(base)[0]
dirname = os.path.dirname(__file__)

outdir = args.o


try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

name = os.path.join(outdir, name_file)

if args.hmm:
    hmm = args.hmm
else:
    hmm = ("/Nancy/mrayko/db/pfam/Pfam-A.hmm")

if args.db:
    blastdb = args.db
else:
    blastdb = ("/Bmo/ncbi_nt_database/nt")


hmmsearch = "/Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch"
prodigal = "/Nancy/mrayko/Libs/Prodigal/prodigal"
cbar = "/Nancy/mrayko/Libs/cBar.1.2/cBar.pl"
#hmm_list_ps01 = "/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/new_table_w_nc/plasmid_hmms_table_ps01.txt"
hmm_list_ps01 = "/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/new_table_w_nc_top_hit/plasmid_hmms_table_ps01.txt"
#list378="/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/378.sorted" 
#list378="/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/378.sorted

# run hmm
#os.system ("prodigal  -p meta -i " + args.f + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log")
#os.system ("hmmsearch  --noali --cut_nc  -o "+name+"_out_pfam --tblout "+name+"_tblout --cpu 10 "+ hmm + " "+name+"_proteins.fa")
#os.system ("tail -n +4 " + name +"_tblout | head -n -10 | awk '{print $1}'| sed 's/_[^_]*$//g'| sort | uniq > " + name +"_plasmid_contigs_names.txt")


# parse hmms
#with open(hmm_list, "r") as infile1:
 #       hmm_list=infile1.readlines()
        
hmm_list = [line.split("\t") for line in open(hmm_list_ps01)]


with open(name + "_tblout", "r") as infile2:
        tblout_pfam=infile2.readlines()



K = 10
if args.K:
  K = int(args.K)

pl_hmms = []
for i in hmm_list:
  if float(i[-1].strip()) >= K:
      pl_hmms.append(i[0])

set378 = set(pl_hmms)


for i in tblout_pfam:
    i=[x for x in i if x]



# Create table backbone.
table =[]

ids=[]
records = list(SeqIO.parse(args.f, "fasta"))
for i in records: 
    ids.append(i.id)  # take all fasta ids and append to table

for i in ids:
  table.append(i.split())



# Collect all hmms for each contig.
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

    item.append(' '.join(hits)) # add plasmid-specific hits
    item.append(' '.join(genes)) # add all hits


    if hits!="-":
        item.append ("hmm+")
    else: 
        item.append ("hmm-")


    if args.c:
   # for i in table:   # print (i)
      if genes!="-":     #i[2]!="-" or i[2][0]!="-":
     #   print (i[2].split(" "))
        print (genes)
        item.append (scikit_multNB(genes))
      else:
         item.append ("Scikit_NBC_Unknown")




if args.b:
    #run blast
 #   os.system ("blastn  -query " + args.f + " -db " + blastdb + " -evalue 0.00001 -outfmt 5 -out "+name+".xml -num_threads 10")
    parser(name+".xml", outdir)
     #    Add blast results
    with open(name + "_span_identity.txt", 'r') as infile:
         span_identity=infile.readlines()
         span_identity = [i.split("\t") for i in span_identity] 
    
    span = {}
    for i in span_identity:
        span[i[0]]=i[1:]

    for item in table:
         #for i in span_identity:
          #   print (item[0])
        # print (i[0]
        if item[0] in span:
            i=span[item[0]]
            item+=[i[3],i[5],i[6].strip(),i[0].strip()]
        else:
           item.append("-")


# Output

with open(name + "_result_table.tsv", "w") as outfile:
    
    if args.b:
        outfile.write('\t'.join(["Contig","Plasmid_hmm", "All_hmm", "Status", "Span", "Identity", "Blast_category", "Best_hit"])+"\n")
    else:
        outfile.write('\t'.join(["Contig","Plasmid_hmm", "Status"])+"\n")

    for i in table:
         outfile.write('\t'.join(i)+"\n")
