import sys


def load_names(filterfile):
    res = set()
    with open(filterfile, "r") as fin:
        for ln in fin.readlines():
            res.add(ln.strip())
    return res

maffile = sys.argv[1]
outprefix = sys.argv[2]
names = set()
if len(sys.argv) > 3:
    filterfile = sys.argv[3]
    names = load_names(filterfile)

res = ""
idealfasta = ""
with open(maffile, "r") as fin:
    s = fin.readline()
    while s != "":
        if s.startswith("a"):
            s1 = fin.readline().strip().split()
            s2 = fin.readline().strip().split()
            name = s2[1]
            if len(sys.argv) <= 3 or name in names:
                start = s1[2]
                l1 = s1[3]
                l2 = s2[3]
                strand = s2[4]
                ref = s1[6].replace("-", "")
                read = s2[6].replace("-", "")
                res += "\t".join([name, start, l1, l2, strand, ref, read]) + "\n"
                idealfasta += ">" + name + "\n" + ref + "\n"
        s = fin.readline()

with open(outprefix + ".tsv", "w") as fout:
    fout.write(res)

with open(outprefix + ".fasta", "w") as fout:
    fout.write(idealfasta)
