import os
import sys
import ntpath

n = int(sys.argv[1])

with open("/Nancy/mrayko/plasmid_prediction_test/downsampling_20perc/chrom_10percen/2260_10fold.txt") as f:
    table = f.read().splitlines()

hmm_list=open(str(n)+"_hmm_list.txt", "w")
for i in table[:n]: 
    print (i.split())
    hmm_list.write(i.split()[0]+"\n")
hmm_list.close()

os.system ("/Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmfetch -f /Nancy/mrayko/db/pfam/Pfam-A.hmm  "+ str(n)+"_hmm_list.txt"+ " > "+str(n)+"_hmm_list.hmm" )
#awk '$6<0.1 {print $1}'  pfam_count_total_fin_sorted10.txt > 2759_10fold_chrom_HMMs.txt
os.system ("/Nancy/mrayko/Libs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmpress "+ str(n)+"_hmm_list.hmm")
