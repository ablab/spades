import os
import sys
import argparse

from operator import itemgetter


def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Median coverage detection by QUAST")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('-f', required = True, help='Input fasta file with contigs')
    parser.add_argument('-r', required = True, help='Reference files')
    parser.add_argument('-o', required = True, help='Output directory')
    return parser.parse_args()


args = parse_args(sys.argv[1:])


hmp_plasmid_refs="/Sid/dantipov/plasmid_data/refs/HMP/"
hmp_contigs="/Sid/dantipov/spades_output/metaplasmid/23.05/HMP/K77/intermediate_contigs.fasta"

jgi_plasmid_refs="/Sid/dantipov/plasmid_data/refs/JGI/"
jgi_contigs="/Sid/dantipov/spades_output/metaplasmid/24.05/JGI/K127/intermediate_contigs.fasta"

synth_plasmid_refs="/Sid/dantipov/plasmid_data/refs/SYNTH/"
synth_contigs="/Sid/dantipov/spades_output/metaplasmid/24.05/SYNTH/K77/intermediate_contigs.fasta"

hmp_chrom_refs="HMP_chrom_refs/"
jgi_chrom_refs="JGI_chrom_refs/"
synth_chom_refs="SYNTH_chrom_refs/"



def run_quast(refs, contigs, out):
    for root, dirs, files in os.walk(refs):  
        for filename in files:
          if (filename.split(".")[-1]) == "fasta" or (filename.split(".")[-1]) == "fa" :
              os.system ("python2 /Nancy/mrayko/Libs/quast-4.5/quast.py -R " + refs + filename + " " + contigs +  " -o " + out  + "/" + filename + " --fast")
    return     


#run_quast(synth_new, synth_contigs, "SYNTH_plasmid_new")
#run_quast(hmp_chrom_refs, hmp_contigs, "HMP_chroms")
#run_quast(jgi_chrom_refs, jgi_contigs, "JGI_chroms")
#run_quast(synth_chrom_refs, synth_contigs, "SYNTH_chroms")



def med_cov (dir):
  for file in next(os.walk(dir))[1]:   
   with open(dir+ "/" + file +"/report.txt") as infile:
     lines = infile.readlines()
     ref_len = lines[18].strip().split(" ")[-1] 


   with open (dir + "/" + file  +"/contigs_reports/nucmer_output/intermediate_contigs.coords.filtered","r") as infile:
   
     hits=infile.readlines()
     hits = [i.strip().split(" ") for i in hits[2:]] 
   long_contigs=[]
   a=0
   for i in hits:
      contig_len=float(i[-1].split("_")[3])
      if abs(float(i[4])-int(i[3])) > (contig_len/2):
         long_contigs.append([contig_len,i[-1].split("_")[-1]])
   sorted_contigs = sorted(long_contigs, key=itemgetter(0), reverse=True)

   for i in sorted_contigs:
     a+=float(i[0])
     if a > float(ref_len)/2:
         print (hits[0][-2],i[1])
         break

  return



def main ():
   run_quast(args.r,args.f,args.o)
   med_cov(args.o)


if __name__ == '__main__':
     main()




