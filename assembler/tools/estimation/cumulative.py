import sys

# check command line arguments
if len(sys.argv) < 3:
	print "Draws cumulative contigs lengths plot"
	print "Usage: ", sys.argv[0], " contigs_file1.fasta contigs_file2.fasta [0 < mul <= 1; default mul = 0.5]"
	exit(0)


if len(sys.argv) >= 4:
	mul = float(sys.argv[3])
else:
	mul = 0.5

# get lengths of contigs from fasta-file #1
lengths1 = []
l = 0
for line in open(sys.argv[1]):
	if (line[0] == '>'):
		if l != 0: # not first sequence in fasta
			lengths1.append(l)
			l = 0
	else:
		l += len(line.strip())
lengths1.append(l)
lengths1.sort(reverse = True)
ln1 = len(lengths1)

# get lengths of contigs from fasta-file #2
lengths2 = []
l = 0
for line in open(sys.argv[2]):
	if (line[0] == '>'):
		if l != 0: # not first sequence in fasta
			lengths2.append(l)
			l = 0
	else:
		l += len(line.strip())
lengths2.append(l)
lengths2.sort(reverse = True)
ln2 = len(lengths2)

# calculate Nx values
vals_percent1 = []
vals_length1 = []
lcur = 0
lind = 0
for l in lengths1:
	lcur += l
	lind += 1
	vals_percent1.append(lind * 100.0 / ln1)
	vals_length1.append(lcur * mul)

vals_percent2 = []
vals_length2 = []
lcur = 0
lind = 0
for l in lengths2:
	lcur += l
	lind += 1
	vals_percent2.append(lind * 100.0 / ln2)
	vals_length2.append(lcur)

# plot vals
import pylab
import matplotlib.ticker
pylab.plot(vals_percent1, vals_length1)
pylab.plot(vals_percent2, vals_length2)
#pylab.yscale('log')
pylab.xlabel('Contigs (percentage)')
pylab.ylabel('Cumulative length')
pylab.title('Cumulative plot')
pylab.grid(True)
ax = pylab.gca()
formatter = matplotlib.ticker.FormatStrFormatter('%.f')
ax.yaxis.set_major_formatter(formatter)
ax.legend((sys.argv[1], sys.argv[2]), loc='lower right')
pylab.savefig('cumulative_plot')
print "Saved to ./cumulative_plot.png"
pylab.show()
