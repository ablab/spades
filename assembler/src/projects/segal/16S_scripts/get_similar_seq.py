import sys
from Bio import SeqIO
import multiprocessing

import editdistance
import edlib

def edist(lst):
    #return editdistance.eval(lst[0], lst[1])
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", additionalEqualities=[('U', 'T')
                                                , ('R', 'A'), ('R', 'G')
                                                , ('Y', 'C'), ('Y', 'T'), ('Y', 'U')
                                                , ('K', 'G'), ('K', 'T'), ('K', 'U')
                                                , ('M', 'A'), ('M', 'C')
                                                , ('S', 'C'), ('S', 'G')
                                                , ('W', 'A'), ('W', 'T'), ('W', 'U')
                                                , ('B', 'C'), ('B', 'G'), ('B', 'T'), ('B', 'U')
                                                , ('D', 'A'), ('D', 'G'), ('D', 'T'), ('D', 'U')
                                                , ('H', 'A'), ('H', 'C'), ('H', 'T'), ('H', 'U')
                                                , ('V', 'A'), ('V', 'C'), ('V', 'G')
                                                , ('N', 'A'), ('N', 'C'), ('N', 'G'), ('N', 'T'), ('N', 'U')] )
    return result["editDistance"]

if len(sys.argv) < 3:
    print "Usage: get_similar_seq.py <fasta wit ideal 16S> <fasta with aligned 16S> "
fasta1 = sys.argv[1]
fasta2 = sys.argv[2]
records = list(SeqIO.parse(fasta1, "fasta"))
best_score = {}
used_aligned = set()
used_ideal = set()
all_aligned = set()
all_ideal = set()
to_remove = set()
for record1 in records:
    all_ideal.add(record1.id)
    to_remove.add(record1)
min_val = 100500
threshold = 100
aligned_rec = {}
cnt = 0
with open(fasta2, "rU") as handle:
    for record2 in SeqIO.parse(handle, "fasta"):
        if cnt % 1000 == 0:
            print cnt, len(to_remove)
        cnt += 1
        all_aligned.add(record2.id)
        seq2 = record2.seq
        aligned_rec[record2.id] = record2
        to_remove_lst = list(to_remove)
        #pool = multiprocessing.Pool(16)
        #vals = pool.map(edist, zip([seq2 for _ in xrange(len(to_remove_lst))], [rec.seq for rec in to_remove_lst]))
        #pool.close()
        #pool.join()
        for i in xrange(len(to_remove_lst)):
            rec = to_remove_lst[i]
            val = edist([seq2, rec.seq])
            if rec.id not in best_score:
                best_score[rec.id] = []
            if val < threshold:
                best_score[rec.id].append(record2)
                to_remove.remove(rec)
                if val < min_val:
                    min_val = val
print "Min ", min_val
res = []
stats = ""
orgs = set()
cnt = 0
print len(best_score.keys())
for it in best_score.keys():
    for rec in best_score[it]:
        used_aligned.add(rec.id)
        used_ideal.add(it)
        #print it, rec.id
        res.append(best_score[it])
        orgs.add(it.split("[")[0])

print "Total ideal = ",len(all_ideal), " aligned = ", len(all_aligned)
print "Mapped ideal = ", len(used_ideal), " aligned = ", len(used_aligned)
#print '\n'.join(all_ideal - used_ideal)
# print "================================================"
# for align in all_aligned - used_aligned:
#     min_val = 100500
#     min_rec = ""
#     for rec in to_remove:
#         val = editdistance.eval(aligned_rec[align].seq, rec.seq)
#         if min_val > val:
#             min_val = val
#             min_rec = rec.id
#     print align, min_rec, min_val