import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_tsv(filename):
    res = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if not ln.startswith("#"):
                query, contig, identity, alignment_length, mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bit_score = ln.strip().split("\t")
                contig_length = int(contig.split('_')[3])
                on_end = 0
                if int(s_start) < 10 or (contig_length - int(s_start)) < 10 or int(s_end) < 10 or (contig_length - int(s_end)) < 10:
                    on_end = 1 
                res.append({"query": query, "contig": contig, "contig_length": contig_length, "identity": float(identity), \
                            "alignment_length": int(alignment_length), "mismatches": int(mismatches), "on_end": on_end, \
                            "gap_opens": int(gap_opens), "q_start": int(q_start), "q_end": int(q_end), \
                            "s_start": int(s_start), "s_end": int(s_end), "evalue": float(evalue), "bit_score": float(bit_score)})
    return res



tsv_table = load_tsv(sys.argv[1])
on_end_cnt = sum([x["on_end"] for x in tsv_table])
print("Number of alignment located on the ends: " + str(int(on_end_cnt*100.0/len(tsv_table))) + "% (" + str(on_end_cnt) + " out of " + str(len(tsv_table)) + ")" )



plt.hist([x["identity"] for x in tsv_table], bins=100)
plt.xlabel("Identity")
plt.savefig("identity.png")

plt.figure()
plt.hist([x["alignment_length"] for x in tsv_table], bins=100)
plt.xlabel("Alignment Length")
plt.savefig("alignment_length.png")

plt.figure()
plt.hist([x["contig_length"] for x in tsv_table], bins=100)
plt.xlabel("Contig Length")
plt.savefig("contig_length.png")