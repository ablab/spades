from Bio.Blast import NCBIXML
import sys
import argparse
from glob import glob
from joblib import Parallel, delayed

# input - BLAST XML of scaffolds.fasta (?)
# output:
# if no alignment at all - nosig
# if good alignment (e<0.001):
#  if "plasmid" in name:
#    <3 matches and close length (0.9) - plasmids
#    others - plasmids_bad
# if "chrom"  in name - chrom
# others - unclass 

def parse_args(args):
    ###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument('--f', help='Path to folder with XML output')
    
    return parser.parse_args()


args = parse_args(sys.argv[1:])


files = glob(str(args.f)+"/*.xml")
print (files)




def parser(f):

    result=open(f,"r")
    nosig = open(f[:-4]+"_no_significant.names", "w")
    chrom = open(f[:-4]+"_chromosome.names", "w")
    plasmids = open(f[:-4]+"_plasmid.names", "w")
    plasmids_bad = open(f[:-4]+"_plasmids_bad.names", "w")
    unclas = open(f[:-4]+"_unclassified.names", "w")
    records= NCBIXML.parse(result)
#item=next(records)
    for item in records:
        print item.query
        pl_len = (int) ((item.query).split('_')[3])
        if len(item.alignments) == 0:
            print "No Significant"
            nosig.write(item.query + '\n')
            continue
        good_aln = []
        for alignment in item.alignments:        
            total_len = 0
            for hsp in alignment.hsps:
                if hsp.expect <0.001:
                    total_len += hsp.align_length
            if (len(alignment.hsps) <= 3 and total_len > 0.9 * pl_len and pl_len > 0.9 * total_len):
                good_aln.append(alignment.title)
        pl = 0
        for seq in good_aln:
            if seq.find('lasmid') != -1:
                plasmids.write(item.query + '\n')
                plasmids.write(seq + '\n')
                pl = 1
                break
        if pl == 0:
            seq = item.alignments[0].title
            if seq.find('lasmid') != -1:
                plasmids_bad.write(item.query + '\n')
                plasmids_bad.write(seq + "\n" + "len " + str (item.alignments[0].length) +" alnments lenghts:" + "\n")
                for hsp in item.alignments[0].hsps:
                    plasmids_bad.write( str (hsp.align_length )+ " ")
                    plasmids_bad.write("\n")
                    pl = 1
            elif (seq.find('hromosome') != -1) or (seq.find("complete genome") != -1):
                chrom.write(item.query + '\n')
                chrom.write(seq + '\n')
            else:
                unclas.write(item.query + '\n')
                unclas.write(item.alignments[0].title + '\n')
    return

Parallel(n_jobs=30)(delayed(parser) (args) for args in files[)
#Parallel(n_jobs=30)(delayed(parser) (args) for args in files)

#                print('score:', hsp.score)
#                print('gaps:', hsp.gaps)
#                print('e value:', hsp.expect)

