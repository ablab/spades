#!/usr/bin/env python

import os, errno
import sys
import argparse
import collections
from math import log
from math import exp
import csv

from parse_blast_xml import parser


def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="HMM-based plasmid verification script")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('-f', required = True, help='Input fasta file')
    parser.add_argument('-o', required = True, help='Output directory')
    parser.add_argument('-b', help='Run BLAST on input contigs', action='store_true')
    parser.add_argument('-db', help='Path to BLAST db')
    parser.add_argument('-hmm', help='Path to Pfam-A HMM database')    
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

ids = []
with open(args.f, "r") as ins:
    for line in ins:
        if line[0]==">":
            ids.append(line[1:].strip())

if args.hmm:
    hmm = args.hmm
else:
    hmm = ("/Nancy/mrayko/db/pfam/Pfam-A.hmm")

if args.db:
    blastdb = args.db
else:
    blastdb = ("/Bmo/ncbi_nt_database/nt")

#hmmsearch="/Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmsearch"
#prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"
#cbar="/Nancy/mrayko/Libs/cBar.1.2/cBar.pl"


# run hmm
print ("Gene prediction...")
os.system ("prodigal -p meta -i " + args.f + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
print ("HMM domains prediction...")
os.system ("hmmsearch  --noali --cut_nc  -o "+name+"_out_pfam --tblout "+name+"_tblout --cpu 10 "+ hmm + " "+name+"_proteins.fa")
print ("Parsing...")
os.system ("tail -n +4 " + name +"_tblout  | head -n -10 | sort -r -k1,1 -k 6,6 | awk '!x[$1]++' > "+name+"_tblout_top_hit" )


tblout_pfam= name + "_tblout" 

def get_table_from_tblout(tblout_pfam):
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()
   
    tblout_pfam = [i.split() for i in tblout_pfam[3:-10]] 
    contigs = collections.OrderedDict()
    for i in range (0, len(tblout_pfam)):
        name = tblout_pfam[i][0].rsplit("_", 1)[0]
        if name not in contigs:
         # if float(tblout_pfam[i][4]) <= 1e-06 :
            contigs[name]=[tblout_pfam[i][2]]
        else:
         # if float(tblout_pfam[i][4]) <= 1e-06:
            contigs[name].append(tblout_pfam[i][2])
           # contigs[name].append(tblout_pfam[i][5])

    out = []
    for key, value in contigs.items():
        out+=[str(key) + " "  +  " ".join(value)]

    return out

tr=os.path.dirname(__file__) + "/plasmid_hmms_table_ps1_top_hit_e06_train.txt"

def naive_bayes(input_list):
    threshold = 0.966999581251
    with open(tr, 'r') as infile:
        table=infile.readlines()
        table = [i.split() for i in table]

# hmm dictionary - for each HMM store plasmid and chromosomal frequency 
    hmm_dict = {}
    for i in table:
      if float(i[5]) >= 10 or float(i[5]) <= 0.1:
        hmm_dict[i[0]] = [float(i[3]),float(i[4])]

# Calculate probabilities for each element of input list
    out_list=[]
    for i in input_list:
        chrom, plasm, chrom_log, plasm_log = 1, 1, 0, 0
        for j in i.split():
          if j in hmm_dict.keys(): 
                plasm=plasm*hmm_dict[j][0]
                plasm_log=plasm_log+log(hmm_dict[j][0])
                chrom=chrom*hmm_dict[j][1]
                chrom_log=chrom_log+log(hmm_dict[j][1])
        if (plasm_log - chrom_log) > threshold: out_list.append(["Plasmid", plasm_log, chrom_log, chrom_log - plasm_log])
        else: out_list.append(["Chromosome",  plasm_log, chrom_log, chrom_log - plasm_log])
     
   
    return out_list 


feature_table = get_table_from_tblout(tblout_pfam) 
feature_table = [i.strip().split(' ', 1) for i in feature_table]

feature_table_names=[]
feature_table_genes=[]
for i in feature_table:
      feature_table_names.append(i[0])
      feature_table_genes.append(i[1])


print ("Classification...")
t=feature_table_genes

k = naive_bayes(t)

names_result={}
for i in range (0,len(k)):
  names_result[feature_table_names[i]] = [k[i][0], feature_table_genes[i]]

final_table=[]
for i in ids:
   if i in names_result:
        final_table.append([i, names_result[i][0],names_result[i][1] ])
   else:
        final_table.append([i, "Chromosome", "--"])


with open(name + '_result_table.csv', 'w') as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(final_table)

print ("Done!")

if args.b:
    #run blast
    os.system ("blastn  -query " + args.f + " -db " + blastdb + " -evalue 0.0001 -outfmt 5 -out "+name+".xml -num_threads 20 -num_alignments 50")
    parser(name+".xml", outdir)
