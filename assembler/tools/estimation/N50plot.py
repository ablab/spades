import sys

# check command line arguments
if len(sys.argv) < 2:
	print "Draws Nx plot (from N01 through N50 to N100)"
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
lengths.sort()
lsum = sum(lengths)

# calculate Nx values
vals_Nx = [0.0]
vals_l = [0]
lcur = 0
for l in lengths:
	lcur += l
	x = lcur * 100.0 / lsum
	vals_Nx.append(vals_Nx[-1] + 1e-10) # eps
	vals_l.append(l)
	vals_Nx.append(x)
	vals_l.append(l)
#print vals_Nx
#print vals_l

# plot vals
import pylab
pylab.plot(vals_Nx, vals_l)
#pylab.yscale('log')
pylab.xlabel('Nx')
pylab.ylabel('Contig length')
pylab.title('Nx plot (N00 to N100)')
pylab.grid(True)
pylab.savefig('Nx_plot')
print "Saved to ./Nx_plot.png"
#pylab.show()
