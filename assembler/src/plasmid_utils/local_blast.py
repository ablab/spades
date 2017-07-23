#! /usr/bin/env python
import sys
import argparse
import os
import subprocess
from joblib import Parallel, delayed
#from familyScreener import familyScreener
#from bitscore import bitscore
#from getFasta import updateFasta
from Bio.Blast import  NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
#file_string = ""    
#x = 1
from glob import glob


def parse_args(args):
    ###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument('--b', help='Path to blast database')
    parser.add_argument('--f', help='Path to folder with SPAdes outputs')   
    return parser.parse_args()


def main():
	args = parse_args(sys.argv[1:])

# get the files

	scaffold_files = glob(str(args.f)+"/*/scaffolds.fasta")
	print (scaffold_files)

# for each scaffolds file - blast online (change later to nr), save as xml
	blastn_args_list = []
	for i in scaffold_files:
		blastn_args_list.append(["blastn","-query", i.strip(), "-db", args.b, "-evalue", "0.001", "-outfmt", "5", "-out",str(i.split("/")[-2])+".xml", "-num_alignments","20"])
		print (blastn_args_list[-1])
		#sequences = open(i).read()
		#sequences = list(SeqIO.parse(i, "fasta"))
		#print (sequences)
		#subprocess.call(blastn_args_list[0])
		#exit(0)
	Parallel(n_jobs=30)(delayed(subprocess.call) (args) for args in blastn_args_list)
		
'''
		blastn_cline = NcbiblastnCommandline(query=i, db=args.b, evalue=0.001, outfmt=5, out=str(i.split("/")[-2])+".xml", num_alignments=20)
		print (blastn_cline)
		blastn_cline
		stdout, stderr = blastn_cline()
'''
		#result_handle = NCBIWWW.qblast("blastn", "nt", sequences, descriptions=10, alignments=10, hitlist_size=1)

#		save_result = open(os.path.join(args.f, str(i.split("/")[1])+".xml"), "w")
		
#		save_result.write(result_handle.read())
#		save_result.close()
#		result_handle.close()



		#outfile.write (i.split("/")[1])
		#finalFile = os.path.join(args.f, str(i.split("/")[1])+".XML")



#for t in sequences: 

#			result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
		#print (t.id + "\n" + t.seq)
		#record = SeqIO.read("m_cold.fasta", format="fasta")
		
		#first_record = next(record_iterator)
		#print(first_record.id)
		#print(first_record.description)


#		second_record = next(record_iterator)
#		print(second_record.id)
#		print(second_record.description)
		#BLAST versus nr

		# name outfile as dir
		


		#outfile.write (i.split("/")[1])
		#finalFile = os.path.join(args.f, str(i.split("/")[1])+".XML")




    #output fasta file!
    #print ("Creating final fasta file")
    #finalFile = os.path.join(fastaPath, "filteredFasta.fa")
    #updateFasta(dfamScreen2, fastafile, finalFile)
    #print ("Files outputted in the following directory: %s" % outputPath)
    #print ("Completed!")






#glob("/path/to/directory/*/")
 #dirs = os.path.join(basePath, "Test", "test.fa")



if __name__ == '__main__':
    main()


#result=open(sys.argv[1],"r")

#fasta_files = [fasta1, fasta2, fasta3, fasta4]
#fasta_files = ["ERR038737_scaffolds"]

#for i in fasta_files:
#    File = open("output"+str(x)+".txt","w")
 #   fasta_string = open(i+".fasta").read() #or make the names fasta1.fasta and just do open(i).read
 #   result_handle = NCBIWWW.qblast("blastn", "nr", fasta_string, hitlist_size=10)


#with open("my_blast.xml", "w") as out_handle:
#		out_handle.write(result_handle.read())

    #blast_records = NCBIXML.parse(result_handle) 
    #or blast_record = NCBIXML.read(result_handle) if you only have one seq in file
    #E_VALUE_THRESH = 0.001
    #for blast_record in blast_records:
     #   for alignment in blast_record.alignments:
      #      for hsp in alignment.hsps:
       #         if hsp.expect < E_VALUE_THRESH:
        ##           file_string += "alignment:",alignment.title+"\n"
          #         file_string += "e-value:",hsp.expect+"\n"
    #x += 1
    #File.write(file_string)

#from Bio.Blast import NCBIWWW
#fasta_string = open("ERR038737_scaffolds.fasta").read()
#result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)
#print (result_handle.read())
