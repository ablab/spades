import os
import sys
import ntpath



#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/chromosome_uniq_hmms.hmm" ### Black list (2208)
#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/chromosome_uniq_genbank_hmms.hmm" ### Deep Black list (508)
#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/non_plasmid_pfamA_hmms.hmm" ### All non-plasmid Pfam (7317)
chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/2260_hmm_list.hmm" ### ### 10 times more in chromosomes
hmmscan=" /Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan" 
prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"


name = os.path.splitext(sys.argv[1])[0]

os.system (prodigal + " -p meta -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
os.system (hmmscan + " --noali --cut_ga  -o "+name+"_out_pfam --tblout "+name+"_tblout --cpu 10 "+ chrom_hmms + " "+name+"_proteins.fa")


with open(name+"_tblout") as f:
    table = f.read().splitlines()

chrom_list=[]
if len(table)>13:
    for i in table[3:-10]:
       if float (i.split()[5]) > 100:
           chrom_list.append(i.split()[2].split('_')[1])

no_singleton_list=[]
for x in chrom_list:
    if chrom_list.count(x)>1:
        no_singleton_list.append(x)

outfile=open (name+"_chromosomal_contigs_names.txt", "w")
for i in list(set(no_singleton_list)): 
   outfile.write(i+"\n")


#os.system ("tail -n +4 "+name+"_tblout | head -n -10 | awk '$5<1e-06 {print $3}'| sed 's/_[^_]*$//g'| sort | uniq > "+name+"_chromosomal_contigs_names.txt")
