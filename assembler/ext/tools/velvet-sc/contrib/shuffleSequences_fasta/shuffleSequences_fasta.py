import sys
from Bio import SeqIO
def interleave(iter1, iter2) :
    while True :
        yield iter1.next()
        yield iter2.next()
f1, f2 = open(sys.argv[1]), open(sys.argv[2])
outfile = open(sys.argv[3], 'w')
format = 'fasta' #or "fastq" or ...
records = interleave(SeqIO.parse(f1, format), SeqIO.parse(f2, format))
count = SeqIO.write(records, outfile, format)
outfile.close()
