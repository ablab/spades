import sys
import yaml

from collections import OrderedDict

from galigner_stats import DataLoader
from galigner_stats import GeneralStatisticsCounter
from galigner_stats import BWAhitsMapper
from galigner_stats import GapsStatistics
from galigner_stats import GapsLengthStatistics
from galigner_stats import GapEditDistanceCounter


class TableValuesCounter:

    def __init__(self, reads, badideal):
        self.total = len(reads) - badideal

    def make_str(self, cnt, total):
            return str(cnt) + " (" + str(cnt*100/total) + "%)"

    def make_str2(self, cnt, cnt_empty, total):
        return str(cnt) + " + " + str(cnt_empty) + " (" + str((cnt + cnt_empty)*100/total) + "%)"

    def get_total_num(self):
        return str(self.total)+ " (100%)"

    def get_mapped_num(self, notmapped):
        return self.make_str(self.total - notmapped, self.total)

    def get_badly_mapped_num(self, path_problems):
        return self.make_str(path_problems, self.total)

    def get_bwagap_problems_num(self, badlyaligned):
        return self.make_str(badlyaligned, self.total)

    def bwa_hits_uncertanity(self, unknown):
        return  self.make_str(unknown, self.total)

    def bwa_hits_failure(self, bwa_problems):
        return self.make_str(bwa_problems, self.total)

    def gaps_failure(self, gaps_cnt_stats):
        return self.make_str2(gaps_cnt_stats["wrong_filled_gaps"], gaps_cnt_stats["empty_gaps"], self.total)

    def ends_failure(self, gaps_cnt_stats):
        return self.make_str2(gaps_cnt_stats["wrong_filled_start"], gaps_cnt_stats["empty_start"], self.total) + " / " +\
                self.make_str2(gaps_cnt_stats["wrong_filled_end"], gaps_cnt_stats["empty_end"], self.total)

    def gap_len(self, gaps_len_stats):
        return str(gaps_len_stats["prefix_len"]) + "/" + str(gaps_len_stats["suffix_len"]) + "/" + str(gaps_len_stats["sum_len"])  

    def get_long_mapping_cnt(self, reads, alignedsubpath):
        cnt = 0
        for r in alignedsubpath.keys():
            p = alignedsubpath[r]
            if (p["mapped_len"])*100/len(reads[r]) > 80:
                cnt += 1
        return self.make_str(cnt,self.total)

    def get_gapped_cnt(self, alignedpaths):
        cnt = 0
        for r in alignedpaths.keys():
            p = alignedpaths[r]
            if p["empty"] > 0:
                cnt += 1
        print cnt
        return self.make_str(cnt,self.total)



def count_table_values(reads, alignedpaths, badideal, notmapped, path_problems, alignedsubpath, badlyaligned,\
                        bwa_problems, gaps_cnt_stats, gaps_len_stats, \
                        gaps_cnt_stats_subpath, gaps_len_stats_subpath, gaps_ed):
    tc = TableValuesCounter(reads, badideal)
    res = OrderedDict([
                 ("Total number of reads", tc.get_total_num()),\
                 ("Mapped with GAligner (#reads)", tc.get_mapped_num(notmapped)),\
                 ("Mapping isn't continiousq (#reads)", tc.get_gapped_cnt(alignedpaths)),\
                 ("Path is not equal to true path (#reads)", tc.get_badly_mapped_num(path_problems)),\
                 ("Path with BWA/Gap problems (#reads)", tc.get_bwagap_problems_num(len(badlyaligned))),\
                 ("Path is wrong. BWA hits uncertainty (#reads)", tc.bwa_hits_uncertanity(gaps_cnt_stats["unknown_num"] - bwa_problems)),\
                 ("Resulting BWA hits failure (#reads)", tc.bwa_hits_failure(bwa_problems)), \
                 ("Gap stage failure (#reads)", tc.gaps_failure(gaps_cnt_stats)), \
                 ("Gap edit distance failure (#Dijkstra_run)", tc.make_str(gaps_ed["failed_ed"], gaps_ed["ed"])), \
                 ("Incorrect prefix/suffix in paths with BWA/Gap fail (#reads)", tc.ends_failure(gaps_cnt_stats)),\
                 ("Median length(in nucs) of skipped prefix/suffix/sum in paths with BWA/Gap fail", tc.gap_len(gaps_len_stats)),\
                 ("Path is subpath of true path (#reads)", tc.get_bwagap_problems_num(len(alignedsubpath))),\
                 ("Path is 80% of true path (#reads)", tc.get_long_mapping_cnt(reads, alignedsubpath)),\
                 ("Incorrect prefix/suffix (#reads) in paths that is subpath", tc.ends_failure(gaps_cnt_stats_subpath)),\
                 ("Median length(in nucs) of skipped prefix/suffix/sum in paths that is subpath", tc.gap_len(gaps_len_stats_subpath))
                 ])

    return res.keys(), res

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

def get_name(path):
    name = path.split("/")[-1]
    res = ""
    if name.startswith("dima_filtering"):
        res = "branch: new_weights; "
    else:
        res = "branch: master; "
    if "_bf_" in name:
        res += "brute_force; "
    else:
        res += "dijkstra; "
    if "_ends" in name:
        res += "restore_ends; "
    if "_ideal_" in name:
        res += "ideal_reads; "
    return name

def save_html(s, fl):
    with open(fl, "w") as fout:
        fout.write(s)


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
    html_name = cfg["print_html"]
    res = {}
    for fl in aligned_files:
        alignedpaths = data_loader.load(fl, "galigner_paths")

        general_stats = GeneralStatisticsCounter(reads, truepaths, alignedpaths)
        badideal = general_stats.cnt_badideal()
        notmapped = general_stats.cnt_notmapped() 
        path_problems = general_stats.cnt_problempaths()
        alignedsubpath, badlyaligned = general_stats.divide_paths()

        print "Total=", len(reads), " ideal=", len(reads) - badideal,  " notmapped=", notmapped
        print "Paths with problems ", path_problems
        print len(alignedsubpath), len(badlyaligned)


        bwahits_mapper = BWAhitsMapper(reads, truepaths, badlyaligned)
        bwa_problems = bwahits_mapper.cnt_problembwa(cfg["print_bwa_hits_failure"])

        gaps_statistics = GapsStatistics(bwahits_mapper, edges, K, cfg["print_gaps_failure"])
        gaps_cnt_stats = gaps_statistics.cnt_wronglyclosedgaps() 

        gaps_ed_counter = GapEditDistanceCounter(bwahits_mapper, edges, K)
        gaps_ed = gaps_ed_counter.count_gap_ed()
        
        gaps_stat_len = GapsLengthStatistics(bwahits_mapper, K)
        gaps_len_stats = gaps_stat_len.cnt_median_alignment_length()

        bwahits_mapper = BWAhitsMapper(reads, truepaths, alignedsubpath)
        gaps_statistics = GapsStatistics(bwahits_mapper, edges, K, cfg["print_gaps_failure"])
        gaps_cnt_stats_subpath = gaps_statistics.cnt_wronglyclosedgaps() 
        gaps_stat_len = GapsLengthStatistics(bwahits_mapper, K)
        gaps_len_stats_subpath = gaps_stat_len.cnt_median_alignment_length()


         
        print "BWA fail ", bwa_problems
        print "Gaps stage problems (in, prefix, suffix, unknown) ", gaps_cnt_stats
        print "Median proportion of length of alignment between two fathest bwa hits to read length", gaps_len_stats

        row_names, res[get_name(fl)] = count_table_values(reads, alignedpaths, badideal, notmapped, path_problems, alignedsubpath, badlyaligned,\
                                                bwa_problems, gaps_cnt_stats, gaps_len_stats, \
                                                gaps_cnt_stats_subpath, gaps_len_stats_subpath,
                                                gaps_ed)

    

    caption_below = ["Read aligned to ref and corresponding ref subseq to graph -- read aligned to reference by BWA MEM and alignment length > 0.8*(read length). After ref subseq mapped to graph -- the result of it is a true path.",\
                     "True path -- path produced by MapRead(), aligning sequence from reference that represents read",\
                     "Resulting BWA hits failure -- BWA hits are on edges that are not in true path or not in correct order",\
                     "Gap stage failure -- gap between two neibouring BWA hits wasn't closed by correct list of edges",\
                     "Gap edit distance failure -- gaps that have shorter ED with ideal path than with found with Dijkstra",\
                     "Prefix -- subpath of edges before first BWA hit edge in the path",\
                     "Suffix -- subpath of edges after last BWA hit edge in the path"]

    table = make_table(res, row_names, caption_below, html_name[:-5])
    save_html(table, "path_statistics_" + html_name)








