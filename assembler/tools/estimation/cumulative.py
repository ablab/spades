import sys
import itertools
import pylab
import matplotlib.ticker
# mailto: vyahhi@gmail.com

# check command line arguments
if len(sys.argv) < 2:
	print "Draws cumulative contigs lengths plot"
	print
	print "Usage:", sys.argv[0], "FASTA1 [MUL1 [FASTA1 [MUL2 ...]"
	print "FASTA is path to .fasta file."
	print "MUL is multiplicator for scaling. Default for last optional mul is 1."
	print
	print "Example: python", sys.argv[0], "../../data/debruijn/we_contigs.fasta 0.5 ../../data/debruijn/velvet_contigs.fa 1"
	print
	exit(0)
if len(sys.argv) % 2 == 0:
	sys.argv.append("1.0")

def get_lengths_from_fastafile(filename):
	lengths = []
	l = 0
	for line in open(filename):
		if (line[0] == '>'):
			if l != 0: # not first sequence in fasta
				lengths.append(l)
				l = 0
		else:
			l += len(line.strip())
	lengths.append(l)
	return lengths

	
lengthses = []
muls = []


for filename, mul in itertools.izip(sys.argv[1::2], sys.argv[2::2]):
	lengthses.append(get_lengths_from_fastafile(filename))
	muls.append(float(mul))

for lengths, mul in itertools.izip(lengthses, muls):
	lengths.sort(reverse = True)
	vals_percent = []
	vals_length = []
	ln = len(lengths)
	lcur = 0
	lind = 0
	for l in lengths:
		lcur += l
		lind += 1
		vals_percent.append(lind * 100.0 / ln)
		vals_length.append(lcur * mul)
	pylab.plot(vals_percent, vals_length)
	
	
#pylab.yscale('log')
pylab.xlabel('Contigs (percentage)')
pylab.ylabel('Cumulative length')
pylab.title('Cumulative plot')
pylab.grid(True)
ax = pylab.gca()
formatter = matplotlib.ticker.FormatStrFormatter('%.f')
ax.yaxis.set_major_formatter(formatter)
ax.legend(sys.argv[1::2], loc='lower right')
pylab.savefig('cumulative_plot')
print "Saved to ./cumulative_plot.png"
pylab.show()
