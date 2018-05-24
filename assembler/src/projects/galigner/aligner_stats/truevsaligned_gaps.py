import yaml
import sys
from truevsaligned_paths import DataLoader
from truevsaligned_paths import GeneralStatisticsCounter
from truevsaligned_paths import BWAhitsMapper
from truevsaligned_paths import GapsStatistics
from truevsaligned_paths import GapsLengthStatistics
import edlib

def edist(lst):
    if len(lst[0]) == 0:
        return len(lst[1])
    if len(lst[1]) == 0:
        return len(lst[0])
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
                                                , ('N', 'C'), ('N', 'C'), ('N', 'G'), ('N', 'T'), ('N', 'U')] )
    return result["editDistance"]


def extract_subpath(subpath, range1, range2, edges, K):
    res = edges[subpath[0]][range1["end"]: len(edges[subpath[0]]) - K]
    if range1["end"] > len(edges[subpath[0]]) - K:
        print "WARNING"
    for ind in xrange(1, len(subpath) - 1):
        res += edges[subpath[ind]][: len(edges[subpath[ind]]) - K]
    res += edges[subpath[len(subpath) - 1]][ : range2["start"]]
    return res


def extract_subseq(seq, range1, range2):
    return seq[range1["end"]: range2["start"]]

if __name__ == "__main__":

    with open(sys.argv[1], 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    data_loader = DataLoader()
    reads = data_loader.load(cfg["reads_path"], "fasta")
    edges = data_loader.load(cfg["edges"], "fasta")
    truepaths = data_loader.load(cfg["truepaths_tsv"], "true_paths")
    aligned_files = data_loader.load(cfg["galigner_out"], "txt")
    print "\n".join(aligned_files)
    K = int(cfg["k"])

    res = {}
    for fl in aligned_files:
        alignedpaths = data_loader.load(fl, "galigner_paths")
        
        bwahits_mapper = BWAhitsMapper(reads, truepaths, alignedpaths)
        for r in alignedpaths.keys():
                if r not in truepaths:
                    continue
                if bwahits_mapper.is_unique(r):
                    true_ind = bwahits_mapper.get_bwa_inds_true(r)
                    aligned_ind = bwahits_mapper.get_bwa_inds_aligned(r)
                    for j in xrange(2, len(true_ind) - 1):
                        if truepaths[r]["path"][true_ind[j - 1]: true_ind[j]] != alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j]]:
                            galigner_seq = extract_subpath(alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j] + 1],\
                                                        alignedpaths[r]["edge_ranges"][aligned_ind[j - 1]], alignedpaths[r]["edge_ranges"][aligned_ind[j]], edges, K)
                            true_seq = extract_subpath(truepaths[r]["path"][true_ind[j - 1]: true_ind[j] + 1],\
                                                        alignedpaths[r]["edge_ranges"][aligned_ind[j - 1]], alignedpaths[r]["edge_ranges"][aligned_ind[j]], edges, K)
                            read_seq = extract_subseq(reads[r], alignedpaths[r]["seq_ranges"][aligned_ind[j - 1]], alignedpaths[r]["seq_ranges"][aligned_ind[j]])
                            if edist([read_seq, true_seq]) < edist([read_seq, galigner_seq]) and len(alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j]+1]) > 2: 
                                print r, edist([read_seq, true_seq]), edist([read_seq, galigner_seq]), edist([galigner_seq, true_seq]),len(read_seq)
                                print truepaths[r]["path"][true_ind[j - 1]: true_ind[j] + 1]
                                print alignedpaths[r]["path"][aligned_ind[j - 1]: aligned_ind[j] + 1]