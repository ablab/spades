import edlib

def load_reads(file):
    res = {}
    name = ""
    with open(file, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith(">"):
                name = ln.strip()[1:]
                res[name] = ""
            else:
                res[name] += ln.strip()
    return res

def load_graph(file):
    res = {}
    name = ""
    with open(file, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith(">"):
                name = ln.strip()[1:]
                res[name] = ""
            else:
                res[name] += ln.strip()
    return res

def get_edgeid(name, rc, edges):
    if not rc:
        return name
    else:
        if str(int(name) + 1) in edges:
            return str(int(name) + 1)
        else:
            return str(int(name) - 1)

def load_path(lst, is_rc, edges):
    for i in xrange(len(lst)):
        edge_id = get_edgeid(lst[i]["name"], is_rc, edges)
        lst[i]["seq"] = edges[edge_id]
        lst[i]["name"] = edge_id
    return lst

def extract_seq(lst):
    res = ""
    for i in xrange(len(lst)):
        it = lst[i]
        #print i, it["name"]
        res += it["seq"][it["e_start"]: it["e_end"] ]
    return res

def extract_ends(lst):
    return lst[0]["r_start"], lst[-1]["r_end"]


def load_alignments(filename):
    res = {}
    fin = open(filename, "r")
    bwa_num = 0
    for ln in fin.readlines():
        cur_read, seq_start, seq_end, rlen, path_dirty, edgelen, score, mapped_seq = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        p_lst = []
        for p in path_dirty.split("]")[:-1]:

            lst = p.split(" ")
            if p.startswith(","):
                edgeid = lst[1]
                lst.pop(0)
            else:
                edgeid = lst[0]
            e_start, e_end, r_start, r_end = int(lst[1].split(",")[0][1:]), int(lst[1].split(",")[1][:-1]), int(lst[2].split(",")[0][1:]), int(lst[2].split(",")[1])
            
            p_lst.append({"name": edgeid, "e_start": e_start, "e_end": e_end, "r_start": r_start, "r_end": r_end})
        res[cur_read] = p_lst
    fin.close()
    return res

def load_idealalignments(filename, reads):
    res = {}
    fin = open(filename, "r")
    bwa_num = 0
    for ln in fin.readlines():
        cur_read, rlen, pathsz, path_dirty, edgelen, edgelen2, mapped_seq = ln.strip().split("\t")
        cur_read = cur_read.split(" ")[0]
        if cur_read in reads:
            p_lst = []
            for p in path_dirty.split("]")[:-1]:
                lst = p.split(" ")
                if p.startswith(","):
                    edgeid = lst[1]
                    lst.pop(0)
                else:
                    edgeid = lst[0]
                e_start, e_end, r_start, r_end = int(lst[1].split(",")[0][1:]), int(lst[1].split(",")[1][:-1]), int(lst[2].split(",")[0][1:]), int(lst[2].split(",")[1])
                p_lst.append({"name": edgeid, "e_start": e_start, "e_end": e_end, "r_start": r_start, "r_end": r_end})
            res[cur_read] = p_lst
    fin.close()
    return res

graph_file = "/Sid/tdvorkina/gralign/C.elegans_synth/graph_idealreads/K55/saves/construction.sqn"
reads_file = "/Sid/tdvorkina/gralign/C.elegans/reads/simulated/sd_10k.fasta" #"/home/tdvorkina/scripts/bwa_test/celegans_ideal2_rc.fasta"


edges = load_graph(graph_file)
reads = load_reads(reads_file)
reads_ideal = load_reads("/home/tdvorkina/scripts/bwa_test/celegans_ideal2_rc.fasta")
alignments = load_alignments("/home/tdvorkina/results/aligner_test/noncanonical_dijkstra_fullpaths_2017-12-27_21-10-32_C.elegans_idealgraph.tsv")
ideal = load_idealalignments("/home/tdvorkina/tmp/algorithmic-biology/assembler/ideal_graph_celegans/construction_noncanonical_ideal_corrected_test2_new.tsv", reads)

reads_paths = {}
reads_rc = {}
nm = "/home/tdvorkina/tmp/algorithmic-biology/assembler/ideal_graph_celegans/construction_noncanonical_corrected_gaps2_test.log"
print nm
with open(nm, "r") as fin:
    for ln in fin.readlines():
        #print ln
        rd, rc, lst = ln.strip().split(";")[0], ln.strip().split(";")[1], ln.strip().split(";")[2:]
        reads_paths[rd] = []
        reads_rc[rd] = False
        if rc == "True":
            reads_rc[rd] = True
        for it in lst:
            filtered = it#.replace(" ", "")
            if ">" in filtered:
                reads_paths[rd].append(filtered.split(">"))
                #print reads_paths[rd][-1][0].split(" ")
                start, end = reads_paths[rd][-1][0].split(" ")[0], reads_paths[rd][-1][0].split(" ")[2]
                reads_paths[rd][-1][0] = [x for x in reads_paths[rd][-1][0].split(" ")[1].split(",") if x != ""]
                reads_paths[rd][-1][1] = [str(int(x)) for x in reads_paths[rd][-1][1].split(",") if x != ""]
                reads_paths[rd][-1].append([int(start), int(end)])
                #print reads_paths[rd][-1]
        #print reads_paths[rd]
cnt = 0
cnt2 = 0
cnt3 = 0
total = 0
for rd in reads_paths.keys():
    if rd not in ideal:
        continue
    la_lst = alignments[rd]
    ideal_lst = ideal[rd]
    read = reads[rd]
    ideal_read = reads_ideal[rd]
    for wp in reads_paths[rd]:
        path1_lst = []
        wasend = False
        s, e = wp[2]
        print s, e
        print wp[0]
        for i in xrange(len(la_lst) - len(wp[0]) + 1):
            j = i
            while j - i < len(wp[0]) and j < len(la_lst) and la_lst[j]["name"] == wp[0][j - i]:
                j += 1
            if j - i == len(wp[0]):
                if s == la_lst[i]["r_start"] and e == la_lst[j - 1]["r_end"]:
                    path1_lst = la_lst[i: j]
                    if j == len(la_lst) - 1:
                        wasend = True
                    break
        path2_lst = []
        print wp[1]
        for i in xrange(len(ideal_lst) - len(wp[1]) + 1):
            j = i
            while j - i < len(wp[1]) and j < len(ideal_lst) and ideal_lst[j]["name"] == wp[1][j - i]:
                j += 1
            if j - i == len(wp[1]): 
                print i, j, ideal_lst[i: j]
                if (s - ideal_lst[i]["r_start"]) * (s - ideal_lst[j - 1]["r_end"]) <= 0 or (s - ideal_lst[i]["r_start"]) * (e - ideal_lst[i]["r_start"]) <= 0 :
                    path2_lst = ideal_lst[i: j]
                    break
        print rd
        print path1_lst
        print path2_lst
        path1_lst = load_path(path1_lst, False, edges)
        path2_lst = load_path(path2_lst, False, edges)

        k = 55
        path1_seq = extract_seq(path1_lst)
        #path2_lst[-1]["e_end"] += k
        path2_seq = extract_seq(path2_lst)
        start1, end1 = extract_ends(path1_lst)
        #print start, end
        if wasend:
           end1 += k 
        subseq1 = read[start1: end1]
        
        ed1 = edlib.align(subseq1, path1_seq, mode="NW" )
        start2, end2 = extract_ends(path2_lst)
        subseq2 = read[start2: end2]
        ed2 = edlib.align(subseq2, path2_seq, mode="NW" )

        if ed1["editDistance"] < ed2["editDistance"]:
            cnt += 1
            # if len(subseq2) < 2000:
            #     print "----", rd, "-----"
            #     print ",".join(wp[0])
            #     print ",".join(wp[1])

            #     print path1_lst
            #     print path2_lst

            #     print start1, end1
            #     print path1_seq
            #     print subseq1
            #     # print ""
            #     print start2, end2
            #     print path2_seq
            #     print subseq2
            #     print "Path1 ed=", ed1["editDistance"]
            #     print "Path2 ed=", ed2["editDistance"]

            #     start, end = extract_ends(path1_lst)
            #     if wasend:
            #        end += k 
            #     subseq1 = ideal_read[start: end]
            #     ed1 = edlib.align(subseq1, path1_seq, mode="NW" )
            #     start, end = extract_ends(path2_lst)
            #     subseq2 = ideal_read[start: end]
            #     ed2 = edlib.align(subseq2, path2_seq, mode="NW" )
            #     print "Path1 ideal ed=", ed1["editDistance"]
            #     print "Path2 ideal ed=", ed2["editDistance"]
            #     print len(read), len(ideal_read)
            #     total += 1
        # if ed2["editDistance"] == 55 or ed2["editDistance"] == 0 :
        #     cnt2 += 1
    # if total > 10:
    #     break

print cnt, total, len(reads_paths.keys())
print cnt2
