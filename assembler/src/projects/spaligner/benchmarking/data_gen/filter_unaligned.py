from Bio import SeqIO
import sys
import os


_, file_extension = os.path.splitext(sys.argv[1])
file_extension = file_extension[1:]
input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), file_extension)
short_seq_iterator = (record for record in input_seq_iterator \
                      if "_aligned" in record.id)

output_handle = open(sys.argv[2], "w")
SeqIO.write(short_seq_iterator, output_handle, file_extension)
output_handle.close()
