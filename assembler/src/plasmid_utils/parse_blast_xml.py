from Bio.Blast import NCBIXML
import sys
import argparse
from glob import glob
#from joblib import Parallel, delayed
import os.path
import pickle


#hsps with larger e-value are discarded
significant_e_value = 0.001

#if query coverage fraction is smaller then alignment is discarded
significant_query_fraction = 0.1

#if best_hit_query_coverage * significant_ratio < other_hit_query_coverage && hits are of different types (i.e plasmid and chromosome then ambiguous)
significant_ratio = 0.9

def parse_args(args):
    ###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument('-i', help='Path to file with XML output')
    parser.add_argument('-o', help='Path to output folder')  
    return parser.parse_args()


###### Return total length of alignment with respect of significant_e_value (0.001) threshold
#
def get_total_len(alignment):
    global significant_e_value
    total_len = 0
    for hsp in alignment.hsps:
        if hsp.expect < significant_e_value:
            total_len += hsp.align_length
    return total_len

def get_query_coverage (alignment):
    global significant_e_value
    total_len = 0
    coords = []
    for hsp in alignment.hsps:
        if hsp.expect < significant_e_value:
            coords.append([hsp.query_start, 1])
            coords.append([hsp.query_end, -1])
    coords.sort()
    if len(coords) > 0:
        layers = coords[0][1]
        for j in range (1, len(coords)):
            if layers > 0:
                total_len  += coords[j][0] - coords[j - 1][0]
            layers += coords[j][1]
           
    return total_len
    

def get_identity(alignment):
    global significant_e_value
    identity = 0
    for hsp in alignment.hsps:
        if hsp.expect < significant_e_value:
            identity += hsp.identities
    return identity

def get_hsp_count(alignment):
    global significant_e_value
    hsp_count = 0
    for hsp in alignment.hsps:
        if hsp.expect < significant_e_value:
            hsp_count += 1
    return hsp_count


def parse_name(seq):
    if (seq.find('hage') != -1 or seq.find('irus') != -1 ):
        return "Virus"
    elif seq.find('lasmid') != -1:
        return "Plasmid"
    elif (seq.find('hromosome') != -1 or seq.find("omplete genome") != -1):
        return "Chromosome"
    else:
        return "Unknown"

def report_file(file, query, alignment):
    file.write(query + '\n')
    file.write(alignment[1] + "\t query coverage: " + str(1.0 * alignment[0]/alignment[6]) + "\t identity: " + str(1.0 * alignment[4]/alignment[3]) + "\t hsp_num: " + str(alignment[5]))
    file.write('\n')
    return



def parser(f, out_dir):
    global significant_query_fraction
    global significant_ratio

    scat_out={}
    xml_file = open(f,"r")
    name = os.path.basename(f)
    name = os.path.join(out_dir, name)
   # print(name)
    nosig = open(name[:-4]+"_no_significant.names", "w")
    chrom = open(name[:-4]+"_chromosome.names", "w")
    plasmids = open(name[:-4]+"_plasmid.names", "w")
    amb_file = open(name[:-4]+"_ambiguous.names", "w")
    records= NCBIXML.parse(xml_file)
    viral = open(name[:-4]+"_viruses.names", "w")
    files = {"Virus": viral, "Plasmid": plasmids, "Chromosome": chrom}
    for item in records:
    #    print (item.query)
	### We are taking query length from contig name provided by SPAdes. Weird, huh?
#        pl_len = (int) ((item.query).split('_')[3])
        pl_len = item.query_length
	###### No alignment - put in non-significant
        if len(item.alignments) == 0:
           # print ("No Significant")
            nosig.write(item.query + '\n\n')
            continue


	###### E > 0.001 
        good_aln = []
        scat_good_aln = []
        #print item.alignments
        alignments = []
        for alignment in item.alignments:
            seq_type = parse_name(alignment.title)
            query_coverage = get_query_coverage(alignment)
            if seq_type == "Unknown" or get_query_coverage(alignment) < significant_query_fraction * pl_len:
                continue
            alignments.append([get_query_coverage(alignment), alignment.title, parse_name(alignment.title), get_total_len(alignment), get_identity(alignment), get_hsp_count(alignment), pl_len])
        if len(alignments)== 0:
           # print ("No Significant")
            nosig.write(item.query + '\n')
            al = item.alignments[0]
            nosig.write(al.title + "\t" + "total query cov " + str(get_query_coverage(al)) + " of " + str(pl_len))
            nosig.write('\n')
        else:
            alignments.sort()
            alignments.reverse()
            type = alignments[0][2]
            best_query_cov = alignments[0][0]
            ambigous = 0
            for i in range (1, len(alignments)):
                if alignments[i][0] < best_query_cov * significant_ratio:
                    break
                if type !=  alignments[i][2]:
                    ambigous = i
                    break
            if ambigous != 0:
               # print ("Ambiguous")
                amb_file.write(item.query + '\n')
                amb_file.write (alignments[0][1] + '\t' + str(best_query_cov) + '\t' + alignments[ambigous][1] + '\t' + str(alignments[ambigous][0]) + '\n')
            else:
               # print (type)
                report_file(files[type], item.query,alignments[0])
    return

def main ():
    parsed_args = parse_args(sys.argv[1:])
    files = glob(str(parsed_args.i)+"/*.xml")
#    print (files)
    parser(parsed_args.i, parsed_args.o)    
#    Parallel(n_jobs=30)(delayed(parser) (infile, parsed_args.o) for infile in files)
#Parallel(n_jobs=30)(delayed(parser) (args) for args in files)

#                print('score:', hsp.score)
#                print('gaps:', hsp.gaps)
#                print('e value:', hsp.expect)


if __name__ == '__main__':
     main()

