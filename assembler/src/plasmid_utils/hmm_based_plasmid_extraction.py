import os, errno
import sys
import ntpath
import argparse
import subprocess
from joblib import Parallel, delayed
from glob import glob
from Bio import SeqIO


name_file = os.path.split(sys.argv[1])[1]

outdir = sys.argv[2]

try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

name=os.path.join(outdir, name_file)

#hmm= "/Nancy/mrayko/db/pfam/Pfam-A.hmm"
hmm="/Nancy/mrayko/db/plasmid_specific_pfam/378_10fold_plasmid_HMMs.hmm"

hmmscan=" /Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan"
prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"


os.system (prodigal + " -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
os.system (hmmscan + " --noali --cut_ga  -o "+name+"_out_pfam --tblout "+name+"_tblout --cpu 10 "+ hmm + " "+name+"_proteins.fa")
os.system ("tail -n +4 " + name +"_tblout | head -n -10 | awk '$5<0.001 {print $3}'| sort | uniq > " + name +"_plasmid_contigs_names.txt")





# get the files
# protein_files = glob(str(args.f)+"/*.fsa")
#protein_files = name + "_proteins.fa"

#handle = open(name + "_proteins.fa", 'r')
#records = list(SeqIO.parse(handle, "fasta"))


# for each protien - run, save as tblout.
#hmmscan_args_list = []
#for record in records:
      
 #     SeqIO.write(record, open(name + "_" + record.id.replace("|", "_")+".fasta", 'w'), "fasta")


#protein_files = glob(name+"_*.fasta")

#print (protein_files) 
# for each protien - run, save as tblout.
#hmmscan_args_list = []
#for i in protein_files:
     # if not os.path.isfile(i.strip()+"_tblout"):
       
 #        hmmscan_args_list.append([hmmscan,"--noali", "--cut_ga", "--cpu", "1"  , "-o", i.strip()+"_pfam", "--tblout", i.strip()+"_pr_tblout", hmm, i.strip()])

#print (hmmscan_args_list[:4])    
#try:
#Parallel(n_jobs=30)(delayed(subprocess.call) (args) for args in hmmscan_args_list)

#except OSError as e: 
 #        if e.errno == 11: 
  #           import time 
   #          print('.') 
    #         time.sleep(0.5) 

#os.system ("cat " + name + " *_pr_tblout   > "+ name +"_total.tblout")
#os.system ("tail -n +4 " + name +"_total_tblout | head -n -10 | awk '$5<0.001 {print $3}'| sed 's/_[^_]*$//g'| sort | uniq > " + name +"_plasmid_contigs_names.txt")

