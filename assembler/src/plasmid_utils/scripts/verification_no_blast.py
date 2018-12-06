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

base = os.path.basename(sys.argv[1])[1]
name_file = os.path.splitext(base)[0]

outdir = sys.argv[2]

try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

name=os.path.join(outdir, name_file)

#hmm= "/Nancy/mrayko/db/pfam/Pfam-A.hmm"
hmm= "/Nancy/mrayko/db/plasmid_specific_pfam/378_10fold_plasmid_HMMs.hmm"
list378="/Nancy/mrayko/PlasmidVerify/plasmid_specific_HMMs/378_hmms.txt" 
blastdb=" /Bmo/ncbi_nt_database/nt"


hmmscan=" /Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan"
prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"
cbar="/Nancy/mrayko/Libs/cBar.1.2/cBar.pl"


# run hmm
#os.system (prodigal + " -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
os.system (hmmsearch + " --noali  -o "+name+"_out_pfam --cut_nc --tblout "+name+"_tblout --cpu 10 "+ hmm + " "+name+"_proteins.fa")
os.system ("tail -n +4 " + name +"_tblout | head -n -10 | awk '{print $3}'| sed 's/_[^_]*$//g'| sort | uniq > " + name +"_plasmid_contigs_names.txt")

# run cbar
#os.system(cbar + " " + sys.argv[1] + " " + name + "_cbar.txt")

# run blast
#os.system ("blastn -query " + sys.argv[1] + " -db " + blastdb + " -evalue 0.00001 -outfmt 5 -out "+name+".xml -num_threads 10")


# parse hmms
with open(list378, "r") as infile1:
        list378=infile1.readlines()

with open(name + "_tblout", "r") as infile2:
        tblout_pfam=infile2.readlines()

list378 = [i.strip() for i in list378] 
for i in tblout_pfam:
    i=[x for x in i if x]


# get list of proteins
pr=[]
proteins = list(SeqIO.parse(name+"_proteins.fa", "fasta"))
for i in proteins: 
   pr.append(i.id)


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
    pr_list=[]
    for j in tblout_pfam[3:-10]:  #check tblout (last 10 - service strings)
        if j[2] and item[0]: # if there is a protein in tblout and SRR in table
            if j[2] not in pr_list: #  if we didn't see this protein yet
                if item[0] in j[2]  and float(j[4]) < 0.01  : # if there's the SRR in protein ID and e-value < 0.01
                    pr_list.append(j[2]) # add this protein to the list for given SRR
                    item[-1].append(str(j[0])) # add to this new field
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

     # add number of total predicted proteins
    counter=0
    for i in pr: 
        if item[0] in i: counter+=1
 
    item.append (str(counter))
     # add number of pfam hits
    item.append(str(len(pr_list)))




# Add HMM results
for i in table:
#    print(i[1])
    if i[1]!="-":
        i.append("hmm+")
    else:
        i.append("hmm-")



# add cbar results 
with open(name + "_cbar.txt", "r") as infile3:
           cbar_res=infile3.readlines()
           
cbar_res = [i.split('\t') for i in cbar_res] 

cbar_list=[]
for i in cbar_res[:10]: 
   if "Plasmid" in i[2]: 
     cbar_list.append(i[0])
     

for i in table:
  if i[0] in cbar_list:
    i.append("cbar+")
  else:
    i.append("cbar-")


# Add classifier

for i in table:   # print (i)
    if i[2]!="-":
        i.append (naive_bayes(i[2].split(" ")))
    else:
         i.append ("NBC_Unknown")



# Add scikit_classifier 
for i in table:   # print (i)
    if i[2]!="-" or i[2][0]!="-":
        print (i[2].split(" "))
        i.append (scikit_multNB(i[2].split(" ")))
    else:
         i.append ("Scikit_NBC_Unknown")




# Output

with open(name + "_result_table.tsv", "w") as outfile:
    for i in table:
         outfile.write('\t'.join(i)+"\n")
