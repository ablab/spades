# Training script for PlasmidVerify
# Input - training set of plasmids and chromosome chunks in fasta format, Pfam-A database
# Required hmmer in path


#1. Import data

import sys, os
import argparse
from operator import itemgetter

threads = str(10)
hmm = "/Nancy/mrayko/db/pfam/Pfam-A.hmm"


def parse_args(args):
###### Command Line Argument Parser
  parser = argparse.ArgumentParser(description="Get feature table for Naive Bayes classifier from training data")
  if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
  parser.add_argument('-pl', required = True, help='Input train plasmid fasta file')
  parser.add_argument('-chr', required = True, help='Input train chromosomes fasta file')
  return parser.parse_args()
 
args = parse_args(sys.argv[1:])

pl_base = os.path.basename(args.pl)
pl_name_file = os.path.splitext(pl_base)[0]

chr_base = os.path.basename(args.chr)
chr_name_file = os.path.splitext(chr_base)[0]


#2. Run hmmscan on train data:


def get_tblout(infile):
   name = os.path.splitext(os.path.basename(infile))[0] 
   print ("Gene prediction...")
   res = os.system ("prodigal -p meta -i " + infile + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
   if res != 0:
     print ("Prodigal run failed")
     exit(1)    
   print ("HMM domains prediction...")
   res = os.system ("hmmsearch  --noali --cut_nc  -o "+name+"_out_pfam --domtblout "+name+"_domtblout --cpu "+ threads + " " + hmm + " "+name+"_proteins.fa")
   if res != 0:
     print ("hmmsearch run failed")
     exit(2)    
   tblout_pfam= name + "_domtblout" 
   return

get_tblout (args.pl)
get_tblout (args.chr)



#3. Train classifier on training data

pfam={}
with open("/Nancy/mrayko/db/pfam/pfam_names.list", 'r') as infile:
    for line in infile:
      pfam[line.strip()]=[0,0]

with open("/Nancy/mrayko/db/pfam/Pfam-A.hmm", 'r') as infile:
    hmms=infile.readlines()
    hmms = [i.strip() for i in hmms] 


#filter by e-value 


with open(pl_name_file + "_domtblout") as f:
    tblout_pl = f.read().splitlines()
    tblout_pl = [i.split() for i in tblout_pl] 
    for i in tblout_pl:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          a = i[3]
          pfam[a][0] += 1

with open(chr_name_file + "_domtblout") as f1:
    tblout_chr = f1.read().splitlines()
    tblout_chr = [i.split() for i in tblout_chr] 
    for i in tblout_chr:
      if i[0][0]!= '#':
        if float(i[6]) <= 1e-06:
          a = i[3]
          pfam[a][1] += 1



# add pseudocounts
pseudo=1.0
table=[]
for key,value in pfam.items():
  if value[0] >= 10 or value[1] >= 10:
    a = (value[0]+pseudo)/len(tblout_pl)
    b = (value[1]+pseudo)/len(tblout_chr)
    table.append([key, value[0], value[1], a, b, a/b])

for i in sorted(table, key=itemgetter(5), reverse = True):
    i = [str(j) for j in i]
    print ('\t'.join(i))
