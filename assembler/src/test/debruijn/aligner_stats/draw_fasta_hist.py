import sys
import pysam
import subprocess
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import SeqIO

def load_fasta(filename):
    record_dict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    return record_dict

mp1 = load_fasta(sys.argv[1])
mp2 = load_fasta(sys.argv[2])

st1 = set(mp1.keys())
st2 = set(mp2.keys())

ln_mapped = []
ln_unmapped = []
for s in st2:
   ln_mapped.append(len(mp1[s].seq))

for s in st1 - st2:
   ln_unmapped.append(len(mp1[s].seq))

plt.hist(ln_mapped, bins=100)
plt.xlim([0, 100000])
plt.savefig("mapped.png")

plt.hist(ln_unmapped, bins=100)
plt.savefig("unmapped.png")