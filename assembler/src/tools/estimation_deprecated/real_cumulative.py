############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Cumulative Contigs Lengths plotter
# mailto: vyahhi@gmail.com

import sys
import itertools
import pylab
import matplotlib.ticker
import fastaparser

# check command line arguments
if len(sys.argv) < 2:
	print "Draws cumulative contigs lengths plot"
	print
	print "Usage:", sys.argv[0], "FASTA1 [MUL1 [FASTA1 [MUL2 ...]"
	print "FASTA is path to .fasta file."
	print "MUL is multiplicator for scaling. Default for last optional mul is 1.0."
	print
	print "Example: python", sys.argv[0], "../../data/debruijn/we_contigs.fasta 0.5 ../../data/debruijn/velvet_contigs.fa 1"
	print
	exit(0)
if len(sys.argv) % 2 == 0: # last default mul = 1.0
	sys.argv.append("1.0")

for filename, mul in itertools.izip(sys.argv[1::2], sys.argv[2::2]):
	# parse
	lengths = fastaparser.get_lengths_from_fastafile(filename)
	lengths.sort(reverse = True)
	# calculate values for the plot
	vals_percent = []
	vals_length = []
	ln = len(lengths)
	lcur = 0
	lind = 0
	for l in lengths:
		lcur += l
		lind += 1
		x = lind  * float(mul)
		vals_percent.append(x)
		y = lcur * float(mul)
		vals_length.append(y)
	# add to plot
	pylab.plot(vals_percent, vals_length)

# customize plot	
#pylab.yscale('log')
pylab.xlabel('Contigs (percentage)')
pylab.ylabel('Cumulative length')
pylab.title('Cumulative plot')
pylab.grid(True)
ax = pylab.gca()
ax.legend(sys.argv[1::2], loc='lower right')
formatter = matplotlib.ticker.FormatStrFormatter('%.f')
ax.yaxis.set_major_formatter(formatter)

# save and show
filename = 'cumulative_plot'
pylab.savefig(filename)
print "Saved to ./" + filename + ".png"
pylab.show()
