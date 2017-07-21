import xlrd
import csv
import argparse
import sys
from glob import glob
import os
import fastaparser
from genericpath import isdir, exists
from os.path import join

parsed_dir = sys.argv[1]
circular_dir = sys.argv[2]
out_dir = sys.argv[3]
os.mkdir(out_dir)



def SRR_to_SRX (SRR):
	os.popen("-db sra -query "+SRR+ " | esummary | xtract -pattern DocumentSummary -element Run@acc Experiment@acc")
	return


def xlsx_to_csv(infile, outfile):
# open the output csv
	with open(infile, 'w') as myCsvfile:
	# define a writer
		wr = csv.writer(myCsvfile, delimiter="\t")

		# open the xlsx file 
		myfile = xlrd.open_workbook(infile)
    	# get a sheet
		mysheet = myfile.sheet_by_index(0)
    # write the rows
		for rownum in range(mysheet.nrows):
			wr.writerow(mysheet.row_values(rownum))
	return



#def create_table():
# given folder - take ERR name, open ERR_plasmid.names
# Contig name, contig seq, len, coverage, Match name,
# open ERR_unclassified.names
# Contig name, contig seq, len, coverage, match to database
#
# input - BLAST XML of scaffolds.fasta (?)
# output:
# if no alignment at all - nosig
# if good alignment (e<0.001):
#  if "plasmid" in name:
#    <3 matches and close length (0.9) - plasmids
#    others - plasmids_bad
# if "chrom" in name - chrom
# others - unclas
# First - take all unclassified, check


files = glob(str(parsed_dir+"/ERR*_unclassified.names"))  

#print (files)
unclassified_contgis=[]


ERRs = []
for file in files:
	#ERRs.append(file[:-16])
	ERRs.append(file.split("/")[-1][:9])


# final_list - list of tuples ERR-name
#print (ERRs)

final_list=[] 

# then open list of circular plasmids and verify.

#print (files)

for i in range (len(ERRs)):
	
	file = files[i]
	#print (file)
	# add 2nd item - list of ids
	#count=0
	with open(file) as f:
		content = f.readlines()
		content = [x.strip() for x in content] 
		for cont in content:
			#print (cont)
		#count+=1
		#if count % 2 != 0: 
			#print (line)  

			final_list.append([ERRs[i][-9:],cont])  # here I need to append contig name


#print (final_list[:10])

for item in final_list:
	#fullname = glob(str(circular_dir+str(item[0])+"_circular.fasta"))
	fullname = join((circular_dir+str(item[0])+"_circular.fasta")) # fasta with circular contigs
#	print (fullname)

	if exists(str(fullname)): 
		contigs = fastaparser.read_fasta(fullname) # open file with circular contigs for given ERR
		for contig in contigs:
			
			if contig[0][1:] == item[1]: 
				print (contig[0][1:])
				print (item[1])
			#	print (join(out_dir + "/" + str(item[0])+"_circular_nonsignificant.fasta"))
				fastaparser.write_fasta_to_file(join(out_dir + "/" + str(item[0])+"_circular_unclassified.fasta"), [contig])

#				print (contig[0])

	#contigs = fastaparser.read_fasta(fullname)


		#else:
		#	entries[i][i+2]=line




	#for string in file:
	#	print (string)
	#run.append (readline(file))
#print (ERRs)
# Contig name, contig seq, len, coverage, Match name,
#plasmids = open(f[:-4]+"_plasmid.names", "w")

#for i in 




#def main():
 #   create_table()

#if __name__ == '__main__':
 #   main()
