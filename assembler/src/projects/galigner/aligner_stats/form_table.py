
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

res = {}
row_names = ["Total number of reads", "Number of aligned reads", "Total reads length (in nucs)", "Aligned length (in nucs)", "Number of non-trivial paths"]
res[align_file] = {"Total number of reads": total_readnum \
                     , "Number of aligned reads": str(aligned_num) + " (" + str("{:.2f}".format(aligned_num*1.0/total_readnum)) + ")" \
                     , "Total reads length (in nucs)": total_readlen\
                     , "Aligned length (in nucs)": str(aligned_len) + " (" + str("{:.2f}".format(aligned_len*1.0/total_readlen)) + ")"
                     , "Number of non-trivial paths": str(path_cnt) + " (" + str("{:.2f}".format(path_cnt*1.0/total_readnum)) + ")"}

table = make_table(res, row_names, caption_below, html_name[:-5])