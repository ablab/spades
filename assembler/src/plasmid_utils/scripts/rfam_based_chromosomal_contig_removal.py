#from Bio.Blast import NCBIXML
import os
import sys
import ntpath
#import argparse
#from glob import glob
#from joblib import Parallel, delayed

# input - BLAST XML of scaffolds.fasta (?)
# output:
# if no alignment at all - nosig
# if good alignment (e<0.001):
#  if "plasmid" in name:
#    <3 matches and close length (0.9) - plasmids
#    others - plasmids_bad
# if "chrom"  in name - chrom
# others - unclass 



#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/chromosome_uniq_hmms.hmm" ### Black list (2208)
#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/chromosome_uniq_genbank_hmms.hmm" ### Deep Black list (508)
#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/non_plasmid_pfamA_hmms.hmm" ### All non-plasmid Pfam (7317)
#chrom_hmms="/Nancy/mrayko/db/chromosomal_proteins_pfam/non_plasmid_pfamA_hmms_7019.hmm ### 7317 without RefSeq plasmid genes
hmm="/Nancy/mrayko/db/rfam/not_in_zymo/total_not_in_zymo.cm"

cmscan="/home/mrayko/.linuxbrew/Cellar/infernal/1.1.2/bin/cmscan"
#hmmscan=" /Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan"
prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"


name = os.path.splitext(ntpath.basename(sys.argv[1]))[0]
#print name

#os.system (prodigal + " -i " + sys.argv[1] + " -a proteins.fa -o genes.fa" )
#os.system (hmmscan + " -o out_pfam --tblout tblout --cpu 10 "+ chrom_hmms + " proteins.fa")
#os.system ("tail -n +4 tblout | head -n -10 | awk '$5<0.001 {print $3}'| sed 's/_[^_]*$//g'| sort | uniq > chromosomal_contigs_names.txt")

#os.system (prodigal + " -i " + sys.argv[1] + " -a "+name+"_proteins.fa -o "+name+"_genes.fa -p meta"  )
os.system (cmscan + " -o "+name+"_out_rfam --tblout "+name+"_rfam_tblout --cpu 10 "+ hmm + " " + sys.argv[1])
os.system ("tail -n +3 "+name+"_rfam_tblout | head -n -10 | awk '$16<1e-06 {print $3}'| sort | uniq > "+name+"_chromosomal_contigs_names.txt")

