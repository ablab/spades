from Bio.Blast import NCBIXML
import sys
import argparse
from glob import glob
from joblib import Parallel, delayed
significant_e_value = 0.001
import os.path

# Input - BLAST XML of scaffolds.fasta (?)
# Output:
# if no alignment at all - non-significant
# if good alignment (e<0.001):
# if "chrom"  in name - chrom
# if "plasmid" in name:
#    <3 matches and close length (0.9) - plasmids
#    others - plasmids_bad
# 
# others - unclass 

def parse_args(args):
    ###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument('-i', help='Path to folder with XML output')
    parser.add_argument('-o', help='Path to output folder')  
    return parser.parse_args()


    ###### Return total length of alignment with respect of significant_e_value (0.001) threshold
def get_total_len(alignment):
    global significant_e_value
    total_len = 0
    for hsp in alignment.hsps:
        if hsp.expect < significant_e_value:
            total_len += hsp.align_length
    return total_len



def parser(f, out_dir):
    xml_file = open(f,"r")
    name = os.path.basename(f)
    name = os.path.join(out_dir, name)
    nosig = open(name[:-4]+"_no_significant.names", "w")
    chrom = open(name[:-4]+"_chromosome.names", "w")
    plasmids = open(name[:-4]+"_plasmid.names", "w")
    plasmids_bad = open(name[:-4]+"_plasmids_bad.names", "w")
    unclas = open(name[:-4]+"_unclassified.names", "w")
    records= NCBIXML.parse(xml_file)

    for item in records:
        print item.query
	### We are taking query length from contig name provided by SPAdes. Weird, huh?
        pl_len = (int) ((item.query).split('_')[3])
	###### No alignment - put in non-significant
        if len(item.alignments) == 0:
            print "No Significant"
            nosig.write(item.query + '\n')
            continue


	###### E > 0.001 
        good_aln = []
        for alignment in item.alignments:        
            total_len = get_total_len(alignment)
            if (len(alignment.hsps) <= 3 and total_len > 0.9 * pl_len and pl_len > 0.9 * total_len):
                good_aln.append(alignment.title)

        pl = 0
        for seq in good_aln:
	    if ((seq.find('hromosome') != -1) or ((seq.find("omplete genome") != -1) and not (seq.find('lasmid') != -1))):
                chrom.write(item.query + '\n')
                chrom.write(seq + '\n')
                pl = 1
                break
	if pl == 0:
		for seq in good_aln:
        	    if seq.find('lasmid') != -1:
                	plasmids.write(item.query + '\n')
	                plasmids.write(seq + '\n')
        	        pl = 1
	                break


        if pl == 0:
            found = False
            for alignment in item.alignments:
                seq = alignment.title
                total_len = get_total_len(alignment)

#alignment is too small to decide
                if total_len < pl_len * 0.1:
                    continue
                if seq.find('lasmid') != -1:
                    plasmids_bad.write(item.query + '\n')
                    plasmids_bad.write(seq + "\n" + "len " + str (alignment.length) +" alnments lenghts:" + "\n")
                    for hsp in alignment.hsps:
                        plasmids_bad.write( str (hsp.align_length )+ " ")
                        plasmids_bad.write("\n")
                    found = True
                    break
                elif (seq.find('hromosome') != -1) or (seq.find("omplete genome") != -1):
                    chrom.write(item.query + '\n')
                    chrom.write(seq + '\n')
                    found = True
                    break
            if not found:
                unclas.write(item.query + '\n')
                unclas.write(item.alignments[0].title + '\n')
    return


parsed_args = parse_args(sys.argv[1:])
files = glob(str(parsed_args.i)+"/*.xml")
print (files)

Parallel(n_jobs=30)(delayed(parser) (infile, parsed_args.o) for infile in files)
#Parallel(n_jobs=30)(delayed(parser) (args) for args in files)

#                print('score:', hsp.score)
#                print('gaps:', hsp.gaps)
#                print('e value:', hsp.expect)

