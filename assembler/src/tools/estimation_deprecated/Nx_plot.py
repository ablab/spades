############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# N00-N50-N100 plotter
# mailto: vyahhi@gmail.com

import sys
import pylab
import matplotlib.ticker
import fastaparser

# check command line arguments
if len(sys.argv) < 2:
	print "Draws Nx plot (from N00 through N50 to N100)"
	print
	print "Usage: python", sys.argv[0], "FASTA1 [FASTA2 [FASTA3 ..."
	print
	print "Example: python", sys.argv[0], "../../data/debruijn/we_contigs.fasta ../../data/debruijn/velvet_contigs.fa"
	exit(0)

for filename in sys.argv[1:]:
	# parse
	lengths = fastaparser.get_lengths_from_fastafile(filename)
	lengths.sort()
	# calculate values for the plot
	vals_Nx = [0.0]
	vals_l = [0]
	lcur = 0
	lsum = sum(lengths)
	for l in lengths:
		lcur += l
		x = lcur * 100.0 / lsum
		vals_Nx.append(vals_Nx[-1] + 1e-10) # eps
		vals_l.append(l)
		vals_Nx.append(x)
		vals_l.append(l)
	# add to plot
	pylab.plot(vals_Nx, vals_l)

# customize plot
pylab.xlabel('Nx')
pylab.ylabel('Contig length')
pylab.title('Nx plot (N00 to N100)')
pylab.grid(True)
ax = pylab.gca()
#ax.legend(["Our Assembler", "Velvet"], loc='lower right')
ax.legend(sys.argv[1:], loc='lower right')
formatter = matplotlib.ticker.FormatStrFormatter('N%.f')
ax.xaxis.set_major_formatter(formatter)

# save and show
filename = 'Nx_plot.pdf'
pylab.savefig(filename)
print "Saved to ./" + filename #+ ".png"
pylab.show()
