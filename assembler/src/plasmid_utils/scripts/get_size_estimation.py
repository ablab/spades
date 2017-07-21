import os
import sys
import time
import pickle
from Bio import SeqIO

max_target = 5
blast_e_value = 0.0001

base = os.path.basename(sys.argv[1])
name_file = os.path.splitext(base)[0]

genomes_len={}
for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
    genomes_len[record.id] = len(record)


res = os.system("prodigal -p meta -i " + sys.argv[1] + " -a " + name_file + "_proteins.fa -o " + name_file + "_prodigal.out 2> "+name_file+"_prodigal.log")
if res != 0:
    print ("Prodigal run failed")
    exit(1)   

os.system("blastp -db /Nancy/mrayko/db/ncbi_viral_proteins/refseq/refseq_viral_proteins -query "+name_file+"_proteins.fa -out "+name_file+"_blast_out.xml  -outfmt 5  -evalue " + str(blast_e_value) + " -max_target_seqs " + str(max_target) + " -num_threads 20")


result_handle = open(name_file+"_blast_out.xml")
from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)


pr_dict={}
genomes = {}
for record in blast_records:  
  genome = record.query.split()[0].rsplit("_", 1)[0]
  if genome not in genomes:
    genomes[genome] = {}

  for alignment in record.alignments:
     if "hypothetical" not in alignment.title.split() and "uncharacterized" not in alignment.title.split():
       pr_name =  record.query.split()[0] 
       pr_dict[alignment.title.split()[1]] = [" ".join(alignment.title.split()[1:]), alignment.title.split("[")[1][:-1]]
       if pr_name in genomes[genome]:
            genomes[genome][pr_name] += [alignment.title.split()[1]]
       else:
            genomes[genome][pr_name] = [alignment.title.split()[1]]


def get_genome_size_for_proteins(proteins):
    #from title get ID 
    from Bio import Entrez
    Entrez.email = "mike.rayko@gmail.com"
    prot_genomes_ids=set()
    for i in proteins:
      handle = Entrez.esearch(db="protein", term=i)
      record = Entrez.read(handle)
      # with protein id get genome id
      genome_record = Entrez.read(Entrez.elink(dbfrom="protein", db = "nuccore", id=record["IdList"][0]))
      time.sleep(1)
      genome_id = (genome_record[0]["LinkSetDb"][0]["Link"][0]["Id"])
      prot_genomes_ids.add(genome_id)
      # using genome id, get genome length
    sizes={}
    for i in prot_genomes_ids:
      handle = Entrez.efetch(db="nuccore", id=i, rettype="docsum")
      record = Entrez.read(handle)
      time.sleep(1)
    #  print (record)
      name = Entrez.read(Entrez.efetch(db = "taxonomy", id=record[0]["TaxId"]))[0]["Lineage"]    
      
      sizes[name]= int(record[0]["Length"])
    return sizes    


def get_genome_size_for_proteins_local(proteins, pr_dict):
    #from title get ID 
    sizes={}
    for i in proteins:
      gen_name = pr_dict[i][1]
      if gen_name in db:
        sizes[str(" ".join(db[gen_name][:3]))]= int(db[gen_name][-1])
    return sizes  

    

import statistics   
import math


with open(os.path.dirname(__file__) + '/db.pkl', 'rb') as f:
    db = pickle.load(f)

for i in genomes:
    print ("="*40)
    print (i)
    sign_proteins = []
    for j in genomes[i]:
      print (j , pr_dict[genomes[i][j][0]][0])
      disp = []
      for k, v in get_genome_size_for_proteins_local(genomes[i][j],pr_dict).items():
            print (k, v)
            disp.append(v)
      if len(disp) > 1:
          CV = statistics.stdev(disp)/statistics.mean(disp)
          print ("CV: ", CV)
          if (CV) < 0.5 :
                sign_proteins.append(statistics.mean(disp)) 
      else:
          if len(disp)>0:
            sign_proteins.append(disp[0])
            
    print ("Sign. genome size estimations: ", sign_proteins)
#    genome_len = int(i.split("_")[3])
    genome_len = genomes_len[i]
    
    if len (sign_proteins) > 1:
      CV = statistics.stdev(sign_proteins)/statistics.mean(sign_proteins)
      mean = statistics.mean(disp)
      stderr = statistics.stdev(sign_proteins)/math.sqrt(len(sign_proteins))
      print ("CV: ", CV)
      if CV < 0.5:
        if  mean*0.8 < genome_len < mean*1.2:
            print ("Full-length")
        else:
            print ("Partial")
      else:
            print ("Variation is too high")
    else:
      if len(sign_proteins) > 0:
        if sign_proteins[0]*0.8 < genome_len < sign_proteins[0]*1.2:
            print ("Full-length")
        else:
            print ("Partial")


