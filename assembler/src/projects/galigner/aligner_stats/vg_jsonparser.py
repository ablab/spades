import json
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
                                                , ('N', 'C'), ('N', 'C'), ('N', 'G'), ('N', 'T'), ('N', 'U')] )
    return result["editDistance"]

def load_edges(gfa_filename):
    res = {}
    with open(gfa_filename, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("S"):
                _, node_id, seq = ln.strip().split("\t")
                res[int(node_id)] = seq
    return res

def load_reads(filename):
    res = {}
    key = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip()
                res[key] = ""
            else:
                res[key] += ln.strip()
    return res

def make_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq):
    if cur_mp["to_length"] == cur_mp["from_length"] and len(cur_mp["sequence"]) == 0:
        graph_seq += ideal_seq[ideal_seq_offset: ideal_seq_offset + int(cur_mp["to_length"])]
    elif cur_mp["to_length"] == cur_mp["from_length"] and len(cur_mp["sequence"]) != 0:
        graph_seq += cur_mp["sequence"]
    elif cur_mp["from_length"] > 0 and cur_mp["to_length"] == 0 and len(cur_mp["sequence"]) != 0:
        graph_seq += cur_mp["sequence"]
    ideal_seq_offset += int(cur_mp["to_length"])
    return ideal_seq_offset, graph_seq

def make_c(s):
    mp = {"A":"T", "T":"A","C":"G", "G":"C"}
    res = "".join([mp[c] for c in list(s)][::-1])
    return res


filename = "/home/tdvorkina/soft/vg/vg_results_ecoli_new2/sd_0001.json"
gfa_filename = "/home/tdvorkina/soft/vg/vg_results_ecoli_new2/ecoli_graph.split.gfa"

reads = load_reads("/Sid/tdvorkina/gralign/E.coli_synth/reads/pacbio/sd_0001.fasta")
edges = load_edges(gfa_filename)

alignments = {}
for k in reads.keys():
    alignments[reads[k]] = {"name": k, "paths": []}

json_file = open(filename, "r")
for ln in json_file.readlines():
    json_data = json.loads(ln)
    ideal_seq = json_data["sequence"]
    if True: #"path" in json_data.keys():
        print ideal_seq
        graph_seq = ""
        graph_seq_all = ""
        ideal_seq_offset = 0
        prev_offset = 0
        prev_node_id = -1
        start_pos = 0
        end_pos = 0
        subpaths = []
        prev_edge = ""
        was_reverse = False
        for item in json_data["path"]["mapping"]:
            if len(item["position"]) != 0:
                #print item["position"]
                node_id, offset = item["position"]["node_id"], 0
                if "offset" in item["position"]:
                    offset = int(item["position"]["offset"]) 
                
                if prev_node_id == -1:
                    start_pos = offset
                elif prev_node_id != node_id:
                    #print prev_node_id, node_id
                    if was_reverse:
                        graph_seq_all += make_c(prev_edge[start_pos:])
                    else:
                        graph_seq_all += prev_edge[start_pos:]
                    subpaths.append(str(prev_node_id))
                    start_pos = 0
		#else:
		    #print offset - prev_offset 
                for edit in item["edit"]:
                    cur_mp = {"to_length": 0, "from_length": 0, "sequence": ""}
                    for k in cur_mp.keys():
                        if k in edit:
                            cur_mp[k] = edit[k]
                    cur_mp["to_length"] = int(cur_mp["to_length"])
                    cur_mp["from_length"] = int(cur_mp["from_length"])
                    ideal_seq_offset, graph_seq = make_edit(cur_mp, ideal_seq_offset, graph_seq, ideal_seq) 
                    offset += int(cur_mp["from_length"])
                prev_offset = offset
                prev_node_id = node_id
                end_pos = offset
                cur_edge = edges[node_id]
                if "is_reverse" in item["position"]:
                    was_reverse = True
                prev_edge = cur_edge
            else:
                graph_seq = ""
                break 
            
        if graph_seq != "":
            if was_reverse:
                graph_seq_all += make_c(prev_edge[start_pos:end_pos])
            else:
                graph_seq_all += prev_edge[start_pos:end_pos]
            subpaths.append(str(prev_node_id))
            
            alignments[ideal_seq]["paths"].append({"len": len(graph_seq), "ed": edist([ideal_seq, graph_seq]), "ed2": edist([ideal_seq, graph_seq_all])})
            #print len(ideal_seq), edist([ideal_seq, graph_seq]), len(json_data["path"]["mapping"])
            # if alignments[ideal_seq]["name"] == "S1_1582":
            #     print ",".join(subpaths)
            #     print graph_seq_all
            #     print ""
            #     print ideal_seq
            #     print ""
            #     print graph_seq
            #     print edist([ideal_seq, graph_seq]), edist([ideal_seq, graph_seq_all]) 
            #     r = alignments[ideal_seq]
            #     print r["name"], len(r["paths"]), len(ideal_seq), len(json_data["path"]["mapping"])
            #print " ".join([str(x["len"]) + ":" + str(x["ed"]) for x in r["paths"]] ) 
            #break


mapped = 0
for r in alignments:
    print r["name"], len(r["paths"])
    if len(r["paths"]) != 0:
        mapped += 1
    print " ".join([str(x["len"]) + ":" + str(x["ed"]) for x in r["paths"]] ) 

print "Total mapped=", mapped
