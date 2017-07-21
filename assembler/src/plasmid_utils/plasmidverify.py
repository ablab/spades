#!/usr/bin/env python

import os, errno
import sys
import argparse
import collections
from math import log
from math import exp
import csv
import operator
import time
import datetime

def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="HMM-based plasmid verification script")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('-f', required = True, help='Input fasta file')
    parser.add_argument('-o', required = True, help='Output directory')
#    parser.add_argument('-b', help='Run BLAST on input contigs', action='store_true')
    parser.add_argument('--db', help='Run BLAST on input contigs with provided database')
    parser.add_argument('--hmm', help='Path to Pfam-A HMM database')    
    parser.add_argument('-t', help='Number of threads')    

    return parser.parse_args()


def get_table_from_tblout(tblout_pfam):
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()
   
    tblout_pfam = [i.split() for i in tblout_pfam[3:-10]]
    for i in tblout_pfam:
        i[13] = float(i[13])

    tblout_pfam.sort(key = operator.itemgetter(0, 13,17), reverse = True)

    top_genes={}
    for i in tblout_pfam:
        if i[0] not in top_genes:
            top_genes[i[0]] = [[i[3],float(i[13]),float(i[17]),float(i[18])]]
        else:
            for j in top_genes[i[0]]:
                start_i, end_i, start_j, end_j = float(i[17]), float(i[18]), float(j[2]), float(j[3])
                 
                if not ((end_i <= start_j) or (start_i >= end_j)):
                    break
                else: 
                    top_genes[i[0]].append([i[3],float(i[13]),start_i,end_i])
                    break


    contigs = collections.OrderedDict()
    for i in top_genes:
        name = i.rsplit("_", 1)[0]
        if name not in contigs:
            contigs[name] = []
            for i in top_genes[i]:
                contigs[name].append(i[0])
        else:
            for i in top_genes[i]:
                contigs[name].append(i[0])

    out = []
    for key, value in contigs.items():
        out+=[str(key) + " "  +  " ".join(value)]

    return out


def naive_bayes(input_list):
    threshold = 0.714327349608 
    tr=os.path.dirname(__file__) + "/plasmid_hmms_table_ps1_top_hit_e06_train.txt"
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
        if (plasm_log - chrom_log) > threshold: out_list.append(["Plasmid", plasm_log, chrom_log, "{0:.2f}".format(plasm_log - chrom_log)])
        else: out_list.append(["Chromosome",  plasm_log, chrom_log, plasm_log - chrom_log])
        
    return out_list 



def main():

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
                ids.append(line.split()[0][1:])
    
    if args.hmm:
        hmm = args.hmm
    else:
        print ("No HMM database provided") 
        exit(1)    
    
    
    if args.db:
        from parse_blast_xml import parser
        blastdb = args.db
    
    if args.t:
        threads = str(args.t)
    else:
        threads = str(20)
    
    
    # run hmm
    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    print ("Gene prediction...")
    res = os.system ("prodigal -p meta -i " + args.f + " -a "+name+"_proteins.fa -o "+name+"_genes.fa 2>"+name+"_prodigal.log" )
    if res != 0:
        print ("Prodigal run failed")
        exit(1)    

    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
    print ("HMM domains prediction...")
    res = os.system ("hmmsearch  --noali --cut_nc  -o "+name+"_out_pfam --domtblout "+name+"_domtblout --cpu "+ threads + " " + hmm + " "+name+"_proteins.fa")
    if res != 0:
        print ("hmmsearch run failed")
        exit(2)    
    print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
   
    print ("Parsing...")
    
    tblout_pfam= name + "_domtblout" 
        
    
    feature_table = get_table_from_tblout(tblout_pfam) 
    feature_table = [i.strip().split(' ', 1) for i in feature_table]
    
    with open(name + '_feature_table.txt', 'w') as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(feature_table)
    
    
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
      names_result[feature_table_names[i]] = [k[i][0],k[i][3], feature_table_genes[i]]
    


    if args.db:
        #run blast
        print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
        print ("Running BLAST...")
    
        os.system ("blastn  -query " + args.f + " -db " + blastdb + " -evalue 0.0001 -outfmt 5 -out "+name+".xml -num_threads "+threads+" -num_alignments 50")
        print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) 
        print ("Parsing BLAST")
        parser(name+".xml", outdir)
        print (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))

    
        #### add blast results
        plasmids= [line.strip().split("\t") for line in open(name + "_plasmid.names")]
    
        plasmids_list={}
        for i in range(0, len(plasmids)-1):
          if len(plasmids[i])==1:
            plasmids_list[plasmids[i][0]] = [float(plasmids[i+1][1].split(":")[1]), float(plasmids[i+1][2].split(":")[1]), plasmids[i+1][0]]
    
        chrom= [line.rstrip().split("\t") for line in open(name + "_chromosome.names")]
        chrom_list={}
        for i in range(0, len(chrom)-1):
          if len(chrom[i])==1:
            chrom_list[chrom[i][0]] = [float(chrom[i+1][1].split(":")[1]), float(chrom[i+1][2].split(":")[1]), chrom[i+1][0]]
    
    
        vir= [line.rstrip().split("\t") for line in open(name + "_viruses.names")]
        vir_list={}
        for i in range(0, len(vir)-1):
          if len(vir[i])==1:
                vir_list[vir[i][0]] = [float(vir[i+1][1].split(":")[1]), float(vir[i+1][2].split(":")[1]), vir[i+1][0]]
    
        nos= [line.rstrip() for line in open(name + "_no_significant.names")]
        nos_list=[]
        for i in nos:
          if len(i.split())==1:
            nos_list.append(i)
        nos_list = [i.strip() for i in nos_list]
        print (len((nos_list)))
    
    
        amb= [line.rstrip() for line in open(name + "_ambiguous.names")]
        amb_list=[]
        for i in amb:
           if i in ids:
            amb_list.append(i)
        amb_list = [i.strip() for i in amb_list] 
    

    
    final_table=[]
    if args.db:
     for i in ids:
       if i in names_result:
            if i in plasmids_list:
              final_table.append([i, names_result[i][0],names_result[i][1], names_result[i][2], "Plasmid", plasmids_list[i][0], plasmids_list[i][1],plasmids_list[i][2]])
            if i in chrom_list:
              final_table.append([i, names_result[i][0],names_result[i][1], names_result[i][2], "Chromosome", chrom_list[i][0], chrom_list[i][1],chrom_list[i][2]])
            if i in vir_list:
              final_table.append([i, names_result[i][0],names_result[i][1], names_result[i][2], "Virus", vir_list[i][0], vir_list[i][1],vir_list[i][2]])
            if i in nos_list:
              print(i)
              final_table.append([i, names_result[i][0],names_result[i][1], names_result[i][2], "Non-significant"])
            if i in amb:
              final_table.append([i, names_result[i][0],names_result[i][1], names_result[i][2], "Ambiguous"])
    
       else:
            final_table.append([i, "Chromosome", "--"])
       print(final_table)
    
    else:
     for i in ids: 
      if i in names_result:
       final_table.append([i, names_result[i][0],names_result[i][1], names_result[i][2]])
      else:
       final_table.append([i, "Chromosome", "--"])
    
    result_file = name + "_result_table.csv"
    with open(result_file, 'w') as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows(final_table)
    

    print ("Done!")
    print ("Verification results can be found in " + os.path.abspath(result_file))
    
    

if __name__ == "__main__":
    main()



