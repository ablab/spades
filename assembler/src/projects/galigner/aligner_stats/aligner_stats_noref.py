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

def is_aligned(al_len, str_len):
    if abs(al_len - str_len) * 1.0 / str_len >= 0.4:
        return False
    else:
        return True

def make_table(results, row_names, caption, name):
    html = """<html><table border="1"><caption>{}</caption><tr><th></th>""".format(name)
    for run_name in sorted(results.keys()):
        html += """<th><div style="width: 200px; height: 50px; overflow: auto">{}</div></th>""".format(run_name)
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
    print "Usage: aligner_stats.py <file with reads .fasta> <file with alignment info ,-separated> <html-file name>"
    exit(-1)

align_files = sys.argv[2]
read_file = sys.argv[1]
html_file = sys.argv[3]
reads_len = load_reads_len(read_file)


print "Reads lens loaded"

total_readlen = 0
total_readnum = 0
for k in reads_len.keys():
    total_readlen += reads_len[k]
    total_readnum += 1

print "Count total len ", total_readlen, total_readnum
res = {}
row_names = ["Total number of reads", "Number of aligned reads", "Total reads length (in nucs)", "Aligned length (in nucs)", "Number of non-trivial paths"]

for align_file in align_files:
    with open(align_file, "r") as fin:
        aligned_len = 0
        aligned_num = 0
        path_cnt = 0
        path_cnt_len = 0
        path_med_len = []
        edge_med_len = []

        aligned_reads = set()
        ln = fin.readline()
        while ln != "":
            [name, start, end, sz, path, path_len, subread, ed] = ln.split("\t")
            name = name.split(" ")[0]
        	if name not in reads_len:
        	    ln = fin.readline()
                continue
            start = int(start)
            end = int(end)
            name = name.split(" ")[0]
            if is_aligned(end-start, reads_len[name]):
                aligned_len += end - start
                aligned_reads.add(name)
                aligned_num += 1
                if len(path_len.split(",")) > 2:
                    path_cnt += 1
                    path_cnt_len += end - start
                    path_med_len.append(len(path_len.split(",")) - 1)
                    edge_med_len.extend([int(x) for x in path_len.split(",")[:-1]])
            ln = fin.readline()

    print "aligned: ", aligned_num, " total: ", total_readnum, " %: ", aligned_num*1.0/total_readnum*100
    print "aligned len: ", aligned_len, " total len: ", total_readlen, " %: ", aligned_len*1.0/total_readlen*100
    print "paths with more than 1 edge: ", path_cnt, " total len: ", path_cnt_len, " %: ", path_cnt*1.0/total_readnum*100
    print "median path len: ", sorted(path_med_len)[len(path_med_len)/2]
    print "median edge len: ", sorted(edge_med_len)[len(edge_med_len)/2]
    print "avg path len: ", sum(path_med_len)/len(path_med_len)
    print "avg edge len: ", sum(edge_med_len)/len(edge_med_len)

    res[align_file] = {"Total number of reads": total_readnum \
                         , "Number of aligned reads": str(aligned_num) + " (" + str("{:.2f}".format(aligned_num*1.0/total_readnum)) + ")" \
                         , "Total reads length (in nucs)": total_readlen\
                         , "Aligned length (in nucs)": str(aligned_len) + " (" + str("{:.2f}".format(aligned_len*1.0/total_readlen)) + ")"
                         , "Number of non-trivial paths": str(path_cnt) + " (" + str("{:.2f}".format(path_cnt*1.0/total_readnum)) + ")"}


caption_below = ["Non-trivial paths -- paths that contain more than 1 edge"]

table = make_table(res, row_names, caption_below, html_name[:-5])
save_html(table, "aligner_stats_no_ref_" + html_name)
