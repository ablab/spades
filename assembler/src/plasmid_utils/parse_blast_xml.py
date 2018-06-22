from Bio.Blast import NCBIXML
import sys
import argparse
from glob import glob
from joblib import Parallel, delayed
significant_e_value = 0.001
import os.path
import pickle

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

def get_identity(alignment):
    global significant_e_value
    identity = 0
    for hsp in alignment.hsps:
        if hsp.expect < significant_e_value:
            identity += hsp.identities
    return identity



def parser(f, out_dir):
    scat_out={}
    xml_file = open(f,"r")
    name = os.path.basename(f)
    name = os.path.join(out_dir, name)
    nosig = open(name[:-4]+"_no_significant.names", "w")
    chrom = open(name[:-4]+"_chromosome.names", "w")
    plasmids = open(name[:-4]+"_plasmid.names", "w")
    plasmids_bad = open(name[:-4]+"_plasmids_bad.names", "w")
    unclas = open(name[:-4]+"_unclassified.names", "w")
    records= NCBIXML.parse(xml_file)
    viral = open(name[:-4]+"_viruses.names", "w")
    for item in records:
    #    print (item.query)
	### We are taking query length from contig name provided by SPAdes. Weird, huh?
#        pl_len = (int) ((item.query).split('_')[3])
        pl_len = item.query_length
	###### No alignment - put in non-significant
        if len(item.alignments) == 0:
            print ("No Significant")
            nosig.write(item.query + '\n')
            continue


	###### E > 0.001 
        good_aln = []
        scat_good_aln = []
        #print item.alignments
        for alignment in item.alignments:        
            total_len = get_total_len(alignment)
            total_identity = get_identity(alignment)            
            if (len(alignment.hsps) <= 3 and total_len > 0.9 * pl_len and pl_len > 0.9 * total_len):
                good_aln.append(alignment.title)
                scat_good_aln.append([str(alignment.title), total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len])                
        pl = 0
        for seq in good_aln:
            if (seq.find('hage') != -1 or seq.find('irus') != -1 ):
                #total_len = get_total_len(item.seq)
                #total_identity = get_identity(item.seq)
                #viral.write(item.query + '\n')
                viral.write(seq + '\n')
                for i in scat_good_aln:
                    if i[0]==seq:
                        scat_out[str(item.query)] = i+["Virus"]
                #scatter_data[item.query]=[total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Virus"]
                pl = 1
                break
        if pl == 0:
            for seq in good_aln:
                if ((seq.find('hromosome') != -1) or ((seq.find("omplete genome") != -1) and not (seq.find('lasmid') != -1))):
                    ##total_len = get_total_len(item.seq)
                    #total_identity = get_identity(item.seq)
                    chrom.write(item.query + '\n')
                    chrom.write(seq + '\n')
                    for i in scat_good_aln:
                      if i[0]==seq:
                        scat_out[str(item.query)] = i+["Chromosome"]
                    #scatter_data[item.query]=[total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Chromosome"]
                    pl = 1
                    break
        if pl == 0:
            for seq in good_aln:
                if seq.find('lasmid') != -1:
                    plasmids.write(item.query + '\n')
                    plasmids.write(seq + '\n')
                    for i in scat_good_aln:
                      if i[0]==seq:
                        scat_out[str(item.query)] = i+["Plasmid"]
                    #scatter_data[item.query]=[total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Plasmid"]
                    pl = 1
                    break


        if pl == 0:
            found = False
            arr = []
            for alignment in item.alignments:
                seq = alignment.title
                total_len = get_total_len(alignment)
                total_identity = get_identity(alignment)            
                arr.append(str(total_len))
#alignment is too small to decide
                if total_len < pl_len * 0.1:
                    continue
                if (seq.find('hage') != -1 or seq.find('irus') != -1 ):
                    viral.write(item.query + '\n')
                    viral.write(seq + '\n')
                    scat_out[str(item.query)] = [alignment.title, total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Virus"]
                    
                    #for i in scatter_data:
                     # if i[1]==seq:
                        #print i 
                      #  scatter_out.append([i]+["Plasmid"])


                    found = True

                    break
                if seq.find('lasmid') != -1:
                    scat_out[str(item.query)] = [alignment.title, total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Plasmid_bad"]
                    plasmids_bad.write(item.query + '\n')
                    plasmids_bad.write(seq + "\n" + "len " + str (alignment.length) +" alnments lenghts:" + "\n")
                    for hsp in alignment.hsps:
                        plasmids_bad.write( str (hsp.align_length )+ " ")
                        plasmids_bad.write("\n")
                    found = True
                    break
                if (seq.find('hromosome') != -1) or (seq.find("omplete genome") != -1):
                    scat_out[str(item.query)] = [alignment.title, total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Chromosome"]
                    chrom.write(item.query + '\n')
                    chrom.write(seq + '\n')                    
                    found = True
                    break
            if not found:


                scat_out[str(item.query)] = [alignment.title, total_len, pl_len, float(total_len)/pl_len, total_identity, float(total_identity)/total_len, "Unclassified"]
                unclas.write(item.query + '\n')
                unclas.write("Best alignment to " + item.alignments[0].title + '\n' + "Total lengths of all alignments:")
                for s in arr:
                    unclas.write(s + " ")
                unclas.write('\n')

    with open (name[:-4]+"_span_identity.txt","w") as span:
        for a,b in scat_out.items():
            span.write(a + "\t" + "\t".join((str(i) for i in b)) + "\n")
#            print a + "\t" + "\t".join((str(i) for i in b))


    #with open('scatter_data_dict.pkl', 'wb') as fid:
     #   pickle.dump(scat_out, fid)    
    #for i in scatter_data:
     #   print str(i).strip('[]')  
#        print str("".join(i[0]+i[2:]))
    return


def main ():
    parsed_args = parse_args(sys.argv[1:])
    files = glob(str(parsed_args.i)+"/*.xml")
#    print (files)

    Parallel(n_jobs=30)(delayed(parser) (infile, parsed_args.o) for infile in files)
#Parallel(n_jobs=30)(delayed(parser) (args) for args in files)

#                print('score:', hsp.score)
#                print('gaps:', hsp.gaps)
#                print('e value:', hsp.expect)


if __name__ == '__main__':
     main()

