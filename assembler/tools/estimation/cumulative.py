import sys

# check command line arguments
if len(sys.argv) < 2:
	print "Draws cumulative contigs lengths plot"
	print "Usage: ", sys.argv[0], " contigs_file.fasta"
	exit(0)
	
# get lengths of contigs from fasta-file
lengths = []
l = 0
for line in open(sys.argv[1]):
	if (line[0] == '>'):
		if l != 0: # not first sequence in fasta
			lengths.append(l)
			l = 0
	else:
		l += len(line.strip())
lengths.append(l)

# prepare lengths
lengths.sort(reverse = True)
ln = len(lengths)

# calculate Nx values
vals_percent = []
vals_length = []
lcur = 0
lind = 0
for l in lengths:
	lcur += l
	lind += 1
	#vals_Nx.append(vals_Nx[-1] + 1e-10) # eps
	#vals_l.append(lcur)
	vals_percent.append(lind * 100.0 / ln)
	vals_length.append(lcur)
#print vals_Nx
#print vals_l

# plot vals
import pylab
import matplotlib.ticker
pylab.plot(vals_percent, vals_length)
#pylab.yscale('log')
pylab.xlabel('Contigs (percentage)')
pylab.ylabel('Cumulative length')
pylab.title('Cumulative plot')
pylab.grid(True)
ax = pylab.gca()
formatter = matplotlib.ticker.FormatStrFormatter('%.f')
ax.yaxis.set_major_formatter(formatter)
pylab.savefig('cumulative_plot')
print "Saved to ./cumulative_plot.png"
pylab.show()
