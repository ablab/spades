import sys

def rc(seq):
	mapping = {"A": "T", "T": "A", "C": "G", "G": "C"}
	res = []
	for c in seq: 
		res.append(mapping[c])
	return "".join(res[::-1])

f = open(sys.argv[1], "r")
reads = {}

for ln in f.readlines():
	if ln.startswith(">"):
		reads[ln.strip()[1:]] = ""
		cur_read = ln.strip()[1:]
	else:
		reads[cur_read] += ln.strip()

f.close()

f = open(sys.argv[1][:-6] + "_rc.fasta" , "w")
with open(sys.argv[2], "r") as f_ind:
    for ln in f_ind.readlines():
        r, isRC = ln.strip().split(" ")
    	f.write(">" + r + "\n")
        if isRC == "True":
	        f.write(rc(reads[r]) + "\n")
        else:
            f.write(reads[r] + "\n")

f.close()

