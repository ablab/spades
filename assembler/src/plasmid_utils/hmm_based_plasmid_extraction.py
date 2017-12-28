#from Bio.Blast import NCBIXML
import os
import sys
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



hmm="/Nancy/mrayko/db/plasmid_specific_pfam/plasmid_hmms_top_24.hmm"
hmmscan=" /Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan"
prodigal="/Nancy/mrayko/Libs/Prodigal/prodigal"



os.system (prodigal + " -p meta -i " + sys.argv[1] + " -a " + sys.argv[1] + "_proteins.fa -o " + sys.argv[1] + "_genes.fa" )
os.system (hmmscan + " -o " + sys.argv[1] + "_out_pfam --tblout " + sys.argv[1]+"_tblout --cpu 10 " + hmm + " " + sys.argv[1] + "_proteins.fa")
os.system ("tail -n +4 " + sys.argv[1]+"_tblout | head -n -10 | awk '$5<0.001 {print $3}'| sed 's/_[^_]*$//g'| sort | uniq > " + sys.argv[1]+"_plasmid_contigs_names.txt")


