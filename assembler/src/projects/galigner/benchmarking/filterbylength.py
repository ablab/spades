from Bio import SeqIO
import sys
import os


input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fastq")
short_seq_iterator = (record for record in input_seq_iterator \
                      if len(record.seq) > int(sys.argv[3]))

output_handle = open(sys.argv[2], "w")
SeqIO.write(short_seq_iterator, output_handle, "fastq")
output_handle.close()