import os
import pickle
import sys
from math import log, exp
from Bio import SeqIO
from parse_blast_xml import parser
blastdb = "/Bmo/ncbi_nt_database/nt"

name = os.path.splitext(os.path.basename(sys.argv[1]))[0]
threads = str(20)

outdir = name+"_blast"

real_len = {}
for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
    real_len[record.id] = len(record)


with open('/Nancy/mrayko/algorithmic-biology/assembler/src/plasmid_utils/scripts/size_estimation_NB/viral_genomes_over5kb.pkl', 'rb') as f:
    genomes = pickle.load(f)

with open('/Nancy/mrayko/algorithmic-biology/assembler/src/plasmid_utils/scripts/size_estimation_NB/viral_genomes_len_over5kb.pkl', 'rb') as f:
    genomes_len = pickle.load(f)

with open('/Nancy/mrayko/algorithmic-biology/assembler/src/plasmid_utils/scripts/size_estimation_NB/proteins_to_genomes_over5kb.pkl', 'rb') as f:
    proteins_to_genomes = pickle.load(f)
    

with open('/Nancy/mrayko/algorithmic-biology/assembler/src/plasmid_utils/scripts/size_estimation_NB/train.pkl', 'rb') as f:
    train = pickle.load(f)




# Gene prediction
print ("Gene prediction...")
#res = os.system ("prodigal -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
#if res != 0:
#   print ("Prodigal run failed")
#   exit(1)   


#Blast
print ("Running BLAST...")
#blastdb = "/Nancy/mrayko/viruses/size_estimation_NB/train_proteins_over5kb"
#os.system ("blastp  -query " + name+"_proteins.fa" + " -db " + blastdb + " -evalue 0.000001 -outfmt 5 -out "+name+".xml -num_threads "+threads+" -num_alignments 5")

#parser(name+".xml", outdir)


#queries =  get_blast_hits(sys.argv[1])



def get_blast_hits(blast_xml):
  #Input - blast xml, output - dict of viruses and corresponding viral hits by proteins
  result_handle = open(blast_xml)
  from Bio.Blast import NCBIXML
  blast_records = NCBIXML.parse(result_handle)
  prot_to_viruses = {}
  for record in blast_records:
      query_id = record.query.split()[0]
      if query_id not in prot_to_viruses:
        prot_to_viruses[query_id] = []
        
      for alignment in record.alignments:
        prot_to_viruses[query_id] += [proteins_to_genomes[alignment.title.split()[1]]]

  virus_to_viruses={}
  for i in prot_to_viruses:
        virus_name = i.split()[0].rsplit("_", 1)[0]
        if virus_name not in virus_to_viruses:
            virus_to_viruses[virus_name] = {}        
        if len(prot_to_viruses[i]) > 0:
          virus_to_viruses[virus_name][i]=set(prot_to_viruses[i])

  return(virus_to_viruses)

queries =  get_blast_hits(name+".xml")

#queries =  get_blast_hits(sys.argv[1])

# Build corresponding probability distribution

eps=0.00001

queries_log_prob = {}

for contig in queries:
 aposteriori = {}
 for k in train:
    aposteriori[k] = []
 for protein in queries[contig]:
   n = len(queries[contig][protein])
   for k in train:
        if k in queries[contig][protein]:
            aposteriori[k] += [1/n - eps]
        else:
            aposteriori[k] += [eps*n/(len(train)-n)] 
    
 p_vir_genes={}          
 for i in aposteriori:
    sum_log = 0
    for prob in aposteriori[i]:
        sum_log += log(prob)
        
    p_vir_genes[i] = sum_log
 
 queries_log_prob[contig] = p_vir_genes

    
def log_to_prob_normalize(input_dict):        
    # Get max
    pr=set()
    for i in input_dict:
        pr.add(input_dict[i])
  #  print (max(pr))
    
    # Subtract max from all logs and exponentiate.
    norm_sum = 0
    for i in input_dict:
        input_dict[i]-=max(pr)
        input_dict[i] = exp(input_dict[i])
        norm_sum += input_dict[i]

    # Normalize
    for i in input_dict:
      input_dict[i] = round(input_dict[i]/norm_sum,10)        
    return input_dict


for i in queries_log_prob:
    queries_log_prob[i] = log_to_prob_normalize(queries_log_prob[i])

    
import operator
for i in queries_log_prob:
    print ("Query: ", i, real_len[i])      
    
    sorted_i = sorted(queries_log_prob[i].items(), key=operator.itemgetter(1), reverse=True)

    out_prob = 0
    for k in sorted_i:
     #   print(i[0])
        if  0.9*real_len[i] < genomes_len[k[0]][0] < 1.1*real_len[i]:
            out_prob += k[1]
     
    if out_prob >= 0.8:
        print ("Full-length", out_prob)
    elif out_prob >= 0.1:
        print ("Variation is too high", out_prob)
    else:
        print ("Partial", out_prob)
 
    
    for k in sorted_i:
        if k[1] > 0.1:
          print (k[0], genomes_len[k[0]][0], genomes_len[k[0]][1], k[1])
    print ("="*40)

