__author__ = 'tanunia'

import sys
import pysam

def load_reads(filename):
    res = {}
    key = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip().split(" ")[0]
                res[key] = ""
            else:
                res[key] += ln.strip()
    return res


def load_reads_len(filename):
    res = {}
    key = ""
    cur_str = ""
    with open(filename, 'r') as infile:
        for ln in infile.readlines():
            if ln.startswith(">"):
                key = ln[1:].strip().split(" ")[0]
                cur_str = ""
                res[key] = 0
            else:
                cur_str += ln.strip()
                res[key] = len(cur_str)

    return res

def load_good_reads(readsbwamem):
    readssam = pysam.AlignmentFile(readsbwamem, "rb")
    mp = {}
    for r in readssam:
        start = r.reference_start
        end = r.reference_end
        name = r.query_name
        if start == -1 or end is None:
            continue
        if abs((end - start) - reads_len[name]) * 1.0 / reads_len[name] < 0.2:
            if name not in mp or abs((end - start) - reads_len[name]) < abs(mp[name]["len"] - reads_len[name]):
                if name not in mp:
                    mp[name] = {}
                mp[name]["len"] = end - start
                mp[name]["start"] = start
                mp[name]["end"] = end
    print "Map loaded #items: ", len(mp)
    names_set = set(mp.keys())
    return names_set

def is_aligned(al_len, str_len):
    if abs(al_len - str_len) * 1.0 / str_len >= 0.4:
        return False
    else:
        return True

def make_table(results, row_names, caption):
    html = """<html><table border="1"><tr><th></th>"""
    for run_name in sorted(results.keys()):
        html += """<th><div style="width: 200px; height: 50px; overflow: auto">{} </div></th>""".format(run_name)
    html += "</tr>"
    for stat in row_names:
        html += "<tr><td>{}</td>".format(stat)
        for run_name in sorted(results.keys()):
            item = results[run_name]
            html += "<td>{}</td>".format(str(item[stat]))
        html += "</tr>"
    html += "</table>"
    html += "<p>{}</p>".format("<br>".join(caption))
    html += "</html>"
    return html

def save_html(s, fl):
    with open(fl, "w") as fout:
        fout.write(s)


if (len(sys.argv) < 4):
    print "Usage: aligner_stats.py <file with .bam> <file with reads .fasta> <file with alignment info> <name>"
    exit(-1)

align_file = sys.argv[3]
read_file = sys.argv[2]
readsbwamem = sys.argv[1]
html_name = sys.argv[4]

reads_len = load_reads_len(read_file)
names_set = load_good_reads(readsbwamem)
print "Reads loaded"

total_readlen = 0
total_readnum = 0
for k in reads_len.keys():
    if k in names_set:
        total_readlen += reads_len[k]
        total_readnum += 1

print "Count total len ", total_readlen, total_readnum

aligned_len = 0
aligned_num = 0
path_cnt = 0
path_cnt_len = 0
path_med_len = []
edge_med_len = []

with open(align_file, "r") as fin:
    ln = fin.readline()
    while ln != "":
        [name, start, end, sz, path, path_len, subread, ed] = ln.split("\t")
        start = int(start)
        end = int(end)
        name = name.split(" ")[0]
        if name in names_set and is_aligned(end-start, reads_len[name]):
            aligned_num += 1
            aligned_len += (end - start)
            if len(path_len.split(",")) > 2:
                path_cnt += 1
                path_cnt_len += sum([int(x) for x in path_len.split(",")[:-1]])
                path_med_len.append(len(path_len.split(",")) - 1)
                edge_med_len.extend([int(x) for x in path_len.split(",")[:-1]])
        ln = fin.readline()
res = {}
row_names = ["Total number of reads", "Number of aligned reads", "Total reads length (in nucs)", "Aligned length (in nucs)", "Number of non-trivial paths"]
res[align_file] = {"Total number of reads": total_readnum \
                     , "Number of aligned reads": str(aligned_num) + " (" + str("{:.2f}".format(aligned_num*1.0/total_readnum)) + ")" \
                     , "Total reads length (in nucs)": total_readlen\
                     , "Aligned length (in nucs)": str(aligned_len) + " (" + str("{:.2f}".format(aligned_len*1.0/total_readlen)) + ")"
                     , "Number of non-trivial paths": str(path_cnt) + " (" + str("{:.2f}".format(path_cnt*1.0/total_readnum)) + ")"}

caption_below = ["Non-trivial paths -- paths that contain more than 1 edge"]

table = make_table(res, row_names, caption_below)
save_html(table, html_name)

print "aligned: ", aligned_num, " total: ", total_readnum, " %: ", aligned_num*1.0/total_readnum*100
print "aligned len: ", aligned_len, " total len: ", total_readlen, " %: ", aligned_len*1.0/total_readlen*100
print "paths with more than 1 edge: ", path_cnt, " total len: ", path_cnt_len, " %: ", path_cnt*1.0/total_readnum*100
print "median path len: ", sorted(path_med_len)[len(path_med_len)/2]
print "median edge len: ", sorted(edge_med_len)[len(edge_med_len)/2]
print "avg path len: ", sum(path_med_len)/len(path_med_len)
print "avg edge len: ", sum(edge_med_len)/len(edge_med_len)
