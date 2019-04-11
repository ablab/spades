#!/usr/bin/python

############################################################################
# Copyright (c) 2015-2017 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import shutil
import getopt
import re
import datetime
import argparse
import subprocess
from traceback import print_exc

sys.path.append('./src/spades_pipeline/')
import process_cfg

#Log class, use it, not print
class Log:

    text = ""

    def log(self, s):
        self.text += s + "\n"    
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()

### Quality assessment ###

# Element of quality assessment map
class MetricEntry:
    config_name = ''
    quast_name = ''
    should_be_higher_than_theshold = True
    is_int = False
    relative_delta = True
    delta_value = 0.1
    simple_parsing = True
    value = 0.0
    assess = True

    def __init__(self, cfg_name, name, higher, is_integer, rela_delta, delta_val, simple_parse = True):
        self.config_name = cfg_name
        self.quast_name = name
        self.should_be_higher_than_theshold = higher
        self.is_int = is_integer
        self.relative_delta = rela_delta
        self.delta_value = delta_val
        if self.should_be_higher_than_theshold:
            self.value = float('inf')
        else:
            self.value = -float('inf')
        self.simple_parsing = simple_parse

    def parse(self, s):
        if self.simple_parsing:
            #E.g. simple number or xxx + yyy part like genes
            return float(s.split()[0])
        else:
            #E.g. reads aligned 100 (100%)
            p = re.search('\((.+)%\)', s).group(1)
            return float(p)
            

# Construct limit map for given metrics ---  list of triplets (config_name, report_name, should_be_higher_that_threshold)
def construct_limit_map(dataset_info, prefix, metric_list, add_all_params = False):
    limit_map = {}
    params = map(lambda x: MetricEntry(prefix + x[0], x[1], x[2], x[3], x[4], x[5]), metric_list)

    if add_all_params or (prefix + 'assess' in dataset_info.__dict__ and dataset_info.__dict__[prefix + 'assess']):
        log.log("Assessing quality results...")
        for p in params:
            if not p.quast_name in limit_map.keys():
                limit_map[p.quast_name] = []

            new_entry = p
            if p.config_name in dataset_info.__dict__:
                #add metric entry that is to be assessed using valued from config
                new_entry.value = float(dataset_info.__dict__[p.config_name])
                limit_map[p.quast_name].append(new_entry)
            elif len(limit_map[p.quast_name]) == 0 or add_all_params:
                #add metric with no value in config, just print it to the log later
                new_entry.assess = False
                limit_map[p.quast_name].append(new_entry)

    return limit_map


# Compare map of metrics with the thesholds
def assess_map(result_map, limit_map):
    res = 0
    for metric in sorted(result_map.keys()):
        if metric in limit_map and len(limit_map[metric]) > 0:
            for entry in limit_map[metric]:
                metric_value = int(result_map[metric]) if entry.is_int else result_map[metric]
                if not entry.assess:
                    log.log(metric + " = " + str(metric_value))
                    continue

                threshold_value = int(entry.value) if entry.is_int else entry.value
                #that metric shouold be higher than threshold (e.g. N50)
                if entry.should_be_higher_than_theshold:
                    if metric_value < threshold_value:
                        log.err(metric + " = " + str(metric_value) + " is less than expected: " + str(threshold_value))
                        res = -1
                    else:
                        log.log(metric + " = " + str(metric_value) + " >= " + str(threshold_value) + " (OK)")
                #that metric shouold be less than threshold (e.g. misassemblies)
                else:
                    if result_map[metric] > entry.value:
                        log.err(metric + " = " + str(metric_value) + " is higher than expected: " + str(threshold_value))
                        res = -1
                    else:
                        log.log(metric + " = " + str(metric_value) + " <= " + str(threshold_value) + " (OK)")
        else:
            log.log(metric + " = " + str(result_map[metric]))

    return res


#parse quast report, returns result_map
def parse_report(report, limit_map):
    columns = []
    values = []
    f = open(report, 'r')
    for l in f:
        row = l.split('\t')
        if len(row) >= 2:
            columns.append(row[0].strip())
            values.append(row[1].strip())
    f.close()

    result_map = {}
    for metric in limit_map.keys():
        if metric in columns:
            result_map[metric] = limit_map[metric][0].parse(values[columns.index(metric)])

    return result_map


# Assess report
def assess_report(report, limit_map, name = ""):
    log.log("Assessing " + name)
    return assess_map(parse_report(report, limit_map), limit_map)


### Reads quality assessment ###

# Construct limit map for reads quality metrics
def construct_reads_limit_map(dataset_info, prefix):
#                                                      metric name, QUAST name, higher than threshold, is int, relative delta, delta value
    return construct_limit_map(dataset_info, prefix, [('min_genome_mapped', "Genome mapped (%)", True, False, False, 0.5),
                                                      ('min_aligned', "Uniquely aligned reads", True, False, False, 0.5)])


# Run reads quality util
def run_reads_assessment(dataset_info, working_dir, output_dir):
    corrected_reads_dataset = os.path.join(output_dir, "corrected/corrected.yaml")
    exit_code = 0

    if not os.path.exists(corrected_reads_dataset):
        log.warn("Corrected reads were not detected in " + corrected_reads_dataset)
        exit_code = 5
    else:
        rq_params = []
        i = 0
        # Setting params
        while i < len(dataset_info.reads_quality_params):
            option = dataset_info.reads_quality_params[i]
            rq_params.append(str(option))
            if i < len(dataset_info.reads_quality_params) - 1 and option == '-r':
                rq_params.append(os.path.join(dataset_info.dataset_path, str(dataset_info.reads_quality_params[i + 1])))
                i += 1
            i += 1

        rq_output_dir = os.path.join(output_dir, "RQ_RESULTS")
        rq_cmd = os.path.join(working_dir, "src/tools/reads_utils/reads_quality.py") + " " + " ".join(rq_params) + " -o " + rq_output_dir + " " + corrected_reads_dataset
        ecode = os.system(rq_cmd)
        if ecode != 0:
            log.warn("Reads quality tool finished abnormally with exit code " + str(ecode))
            exit_code = 6
        else:
            # Assessing reads quality report
            limit_map = construct_reads_limit_map(dataset_info, "")
            if assess_report(os.path.join(rq_output_dir, "report.tsv"), limit_map) != 0:
                exit_code = 7

    return exit_code


### QUAST ###

# Run QUAST for a set of contigs
def run_quast(dataset_info, contigs, quast_output_dir, opts):
    if not reduce(lambda x, y: os.path.exists(y) and x, contigs, True):
        log.err("No contigs were found")
        return 8

    cmd = "quast.py"
    if dataset_info.mode == "meta":
        cmd = "metaquast.py"
    elif dataset_info.mode == "rna":
        cmd = "rnaQUAST.py"
    #Preparing params
    path_options = ['-R', '-G', '-O']
    if dataset_info.mode == "rna":
        path_options = ['-r', '--reference', '--gmap_index', '--gene_db', '--gtf']

    quast_params = [opts]
    if 'quast_params' in dataset_info.__dict__ and dataset_info.quast_params:
        i = 0
        while i < len(dataset_info.quast_params):
            option = dataset_info.quast_params[i]
            quast_params.append(str(option))
            if i < len(dataset_info.quast_params) - 1 and option in  path_options:
                quast_params.append(os.path.join(dataset_info.dataset_path, str(dataset_info.quast_params[i + 1])))
                i += 1
            i += 1

    log.log('Running ' + cmd + ' on ' + ','.join(contigs))
    quast_home = os.environ['QUAST_HOME'] if "QUAST_HOME" in os.environ and os.environ['QUAST_HOME'] != "" else dataset_info.quast_dir
    quast_cmd = os.path.join(quast_home, cmd) + " " + " ".join(quast_params) + " "

    ctg_option = ""
    if dataset_info.mode == "rna":
        ctg_option = " -c "

    quast_cmd += " -o " + quast_output_dir + " " + ctg_option + " ".join(contigs) + " > /dev/null"
    log.log('Executing ' + quast_cmd)
    ecode = os.system(quast_cmd)
    if ecode != 0:
        log.err(cmd + " finished abnormally with exit code " + str(ecode))
        return 9

    return 0


# Construct limit map for QUAST metrics
def construct_quast_limit_map(dataset_info, prefix, add_all_params = False):
#                                       metric name,  QUAST name, higher than threshold, is int, relative delta, delta value
    return construct_limit_map(dataset_info, prefix, [
                                        ('min_n50',         "N50",                      True,   True, True, 0.1),
                                        ('max_n50',         "N50",                      False,  True, True, 0.1),
                                        ('min_ng50',        "NG50",                     True,   True, True, 0.1) ,
                                        ('max_ng50',        "NG50",                     False,  True, True, 0.1),
                                        ('min_lg50',        "LG50",                     True,   True, True, 0.05) ,
                                        ('max_lg50',        "LG50",                     False,  True, True, 0.05),
                                        ('max_mis',         "# misassemblies",          False,  True, False, 1),
                                        ('min_mis',         "# misassemblies",          True,   True, False, 1),
                                        ('min_genome_mapped', "Genome fraction (%)",    True,   False, False, 0.5),
                                        ('min_genes',       "# genes",                  True,   True, True, 0.01),
                                        ('max_indels',      "# indels per 100 kbp",     False,  False, False, 1),
                                        ('max_subs',        "# mismatches per 100 kbp", False,  False, False, 1),
                                        ('max_localmis',    "# local misassemblies",    False,  True, False, 1),
                                        ('max_ns',          "# N's per 100 kbp",        False,  False, True, 0.05),
                                        ('max_dr',          "Duplication ratio",        False,  False, False, 0.03)], add_all_params)

# Construct limit map for rnaQUAST metrics
def construct_rnaquast_limit_map(dataset_info, prefix, add_all_params = False):
#                                       metric name, QUAST name, higher than threshold, is int, relative delta, delta value
    return construct_limit_map(dataset_info, prefix, [
                                        ('min_transcripts',     "Transcripts",          True,   True, True, 0.1),
                                        ('min_transcripts_500', "Transcripts > 500 bp", True,   True, True, 0.05),
                                        ('min_aligned',         "Aligned",              True,   True, True, 0.05),
                                        ('max_unaligned',       "Unaligned",            False,  True, True, 0.05),
                                        ('min_db_cov',          "Database coverage",    True,   False, True, 0.1),
                                        ('max_db_cov',          "Database coverage",    False,  False, True, 0.1),
                                        ('min_50_genes',        "50%-assembled genes",  True,   True, True, 0.05),
                                        ('min_95_genes',        "95%-assembled genes",  True,   True, True, 0.05),
                                        ('max_95_genes',        "95%-assembled genes",  False,  True, True, 0.05),
                                        ('min_95_cov_genes',    "95%-covered genes",    True,   True, True, 0.05),
                                        ('max_95_cov_genes',    "95%-covered genes",    False,   True, True, 0.05),
                                        ('min_50_iso',          "50%-assembled isoforms", True,  True, True, 0.05),
                                        ('min_95_iso',          "95%-assembled isoforms", True,  True, True, 0.05),
                                        ('max_95_iso',          "95%-assembled isoforms", False, True, True, 0.05),
                                        ('max_mis',             "Misassemblies",        False,   True, True, 0.05)], add_all_params)


# Run QUAST and assess its report for a single contig file
def quast_run_and_assess(dataset_info, fn, output_dir, name, prefix, special_exit_code, opts):
    if os.path.exists(fn):
        log.log("Processing " + fn)
        qcode = run_quast(dataset_info, [fn], output_dir, opts)
        if qcode != 0:
            log.err("Failed to run QUAST!")
            return qcode

        limit_map = None
        if dataset_info.mode == "rna":
            limit_map = construct_rnaquast_limit_map(dataset_info, prefix)
        else:
            limit_map = construct_quast_limit_map(dataset_info, prefix)

        report_path = output_dir
        if dataset_info.mode == "meta":
            report_path = os.path.join(report_path, "combined_reference")

        if dataset_info.mode == "rna":
            report_path = os.path.join(report_path, "short_report.tsv")
        else:
            report_path = os.path.join(report_path, "report.tsv")

        if assess_report(report_path, limit_map, name) != 0:
            return special_exit_code

        return 0
    else:
        log.err("File not found " + fn)
        return 8


# Run QUAST on all contigs that need to be evaluated
def quast_analysis(contigs, dataset_info, folder):
    exit_code = 0
    for name, file_name, prefix, opts, ext in contigs:
        if ext != "fasta":
            continue
        log.log("======= " + name.upper() + " SUMMARY =======")
        if prefix != "":
            prefix = prefix + "_"
        qcode = quast_run_and_assess(dataset_info, os.path.join(folder, file_name + "." + ext),
                        os.path.join(folder, "QUAST_RESULTS_" +name.upper()), name, prefix, 10, opts)
        if qcode != 0:
            log.err("Analysis for " + name + " did not pass, exit code " + str(qcode))
            exit_code = qcode

    return exit_code


### Compare misassemblies 

def parse_alignment(s): 
    m = s.strip()
    if not m.startswith("Real Alignment"):
        return None
    
    coords = m.split(':')[1]
    genome_c = coords.split('|')[0].strip().split(' ')
    genome_coords = (int(genome_c[0]), int(genome_c[1]))
    contig_c = coords.split('|')[1].strip().split(' ')
    contig_coords = (int(contig_c[0]), int(contig_c[1]))
    
    return (genome_coords, contig_coords)


def parse_misassembly(m1, m2):
    m1_coords = parse_alignment(m1)
    m2_coords = parse_alignment(m2)

    if not m1_coords or not m2_coords:
        return None

    pos1 = 0
    if (m1_coords[1][0] > m1_coords[1][1]):
        pos1 = m1_coords[0][0]
    else:
        pos1 = m1_coords[0][1]

    pos2 = 0
    if (m2_coords[1][0] < m2_coords[1][1]):
        pos2 = m2_coords[0][0]
    else:
        pos2 = m2_coords[0][1]

    if pos1 < pos2:
        return (pos1, pos2)
    else:
        return (pos2, pos1)
  

def find_mis_positions(contig_report):
    infile = open(contig_report)
    #skipping prologue
    line = infile.readline()
    while not line.startswith("Analyzing contigs..."):
        line = infile.readline()

    # main part of quast output
    prev_line = ""
    coords = {}
    currend_id = ""
    while not line.startswith("Analyzing coverage..."):
        line = infile.readline()
        if line.startswith("CONTIG:"):
            currend_id = line.split()[1].strip()
        elif (line.find("Extensive misassembly") != -1):
            line2 = infile.readline()
            if (line2.find("Real Alignment") == -1):
                line2 = infile.readline()
            mis = parse_misassembly(prev_line, line2)
            if not mis:
                log.warn("Failed to parse misassembly")
            elif mis not in coords:
                coords[mis] = [currend_id]
            else:
                coords[mis].append(currend_id)
            prev_line = line2
            continue
        prev_line = line
               
    infile.close()
    return coords


# Compare two contig reports from QUAST
def cmp_misassemblies(quast_output_dir, old_ctgs, new_ctgs):

    log.log("Comparing misassemblies in " + old_ctgs + " and " + new_ctgs)

    old_contigs_report = os.path.join(quast_output_dir, "contigs_reports/contigs_report_" + old_ctgs + ".stdout")
    new_contigs_report = os.path.join(quast_output_dir, "contigs_reports/contigs_report_" + new_ctgs + ".stdout")

    if not os.path.exists(old_contigs_report):
        log.warn("Old contigs report was not found in " + old_contigs_report)
        return True
    if not os.path.exists(new_contigs_report):
        log.warn("New contigs report was not found in " + new_contigs_report)
        return True

    old_pos = find_mis_positions(old_contigs_report)
    new_pos = find_mis_positions(new_contigs_report)

    for k, v in old_pos.items():
        if k in new_pos:
            v2 = new_pos[k]
            if len(v) == len(v2):
                old_pos[k] = []
                new_pos[k] = []

    old_mis = 0
    for k, v in old_pos.items(): 
        old_mis += len(v)
    new_mis = 0
    for k, v in new_pos.items():   
        new_mis += len(v)

    if new_mis == 0 and old_mis == 0:
        log.log("Misassemlies are all the same")
        return True
    else:
        log.log("Detected misassembly changes")
        if old_mis == 1:
            log.log(str(old_mis) + " old misassemly is gone:")
        else:
            log.log(str(old_mis) + " old misassemlies are gone:")
        for k, v in old_pos.items():
            if len(v) > 0:
                log.log("   Misassembly genomic positions " + str(k) + " found in " + str(len(v)) + " contig(s)")
                for ctg in v:
                    log.log("      " + ctg)

        if new_mis == 1:
            log.log(str(new_mis) + " new misassembly was detected:")
        else:
            log.log(str(new_mis) + " new misassemblies were detected:")
        for k, v in new_pos.items():
            if len(v) > 0:
                log.log("   Misassembly genomic positions " + str(k) + " found in " + str(len(v)) + " contig(s)")
                for ctg in v:
                    log.log("      " + ctg)
        return False

        
# Run QUAST and compare misassemblies for given pair of contigs files
def compare_misassemblies(contigs, dataset_info, contig_storage_dir, output_dir):
    exit_code = 0
    rewrite_latest = True
    enable_comparison = (dataset_info.mode in ("standard", "tru", "meta", "plasmid"))

    if enable_comparison and contig_storage_dir != '' and 'quast_params' in dataset_info.__dict__ and dataset_info.quast_params and '-R' in dataset_info.quast_params:
        for name, file_name, prefix, opts, ext in contigs:
            if ext != "fasta":
                continue
            fn = os.path.join(output_dir, file_name + "." + ext) 
            if not os.path.exists(fn):
                log.err('File for comparison is not found: ' + fn)
                exit_code = 8

            latest_ctg = os.path.join(contig_storage_dir, "latest_" + name + "." + ext)
            if os.path.exists(latest_ctg):
                log.log("======= " + name.upper() + " COMPARISON =======")
                quast_output_dir = os.path.join(output_dir, "QUAST_RESULTS_CMP_" + name.upper())
                qcode = run_quast(dataset_info, [fn, latest_ctg], quast_output_dir, opts)


                if dataset_info.mode == "meta":
                    quast_output_dir = os.path.join(quast_output_dir, "combined_reference")

                if not cmp_misassemblies(quast_output_dir, "latest_" + name, file_name):
                    rewrite_latest = False
            else:
                log.warn('Failed to find ' + latest_ctg + ', nothing to compare misassemblies with')

    return exit_code, rewrite_latest


### Support functions ###

# Load config
def load_info(dataset_path):
    if not os.path.exists(dataset_path):
        log.err(dataset_path + " can not be found")
        sys.exit(-2)

    if os.path.isdir(dataset_path):
        dataset_path = os.path.join(dataset_path, "dataset.info")

    info = process_cfg.load_config_from_file(dataset_path)
    info.__dict__["dataset_path"] = os.path.split(dataset_path)[0]

    if "mode" not in info.__dict__:
        if 'truseq' in info.__dict__ and info.truseq:
            info.__dict__["mode"] = "tru"
        elif 'meta' in info.__dict__ and info.meta:
            info.__dict__["mode"] = "meta"
        elif 'rna' in info.__dict__ and info.rna:
            info.__dict__["mode"] = "rna"
        elif 'plasmid' in info.__dict__ and info.plasmid:
            info.__dict__["mode"] = "plasmid"
        else:
            info.__dict__["mode"] = "standard"
    return info


# Get contig list to be assessed for this run
def get_contigs_list(args, dataset_info, before_rr = False):
    contigs = [("contigs", "contigs", "", "", "fasta")]
    contigs.append(("scaffolds", "scaffolds", "sc", " -s ", "fasta"))
    contigs.append(("contigs_paths", "contigs", "", "", "paths"))
    contigs.append(("scaffolds_paths", "scaffolds", "", "", "paths"))

    if dataset_info.mode == "tru":
        contigs = [("contigs", "truseq_long_reads", "", "", "fasta")]
    if dataset_info.mode == "rna":
        contigs = [("transcripts", "transcripts", "", "", "fasta")]
    if dataset_info.mode == "meta":
        contigs.append(("preliminary", "first_pe_contigs", "prelim", "", "fasta"))
    if before_rr and dataset_info.mode in ("standard", "meta"):
        contigs.append(("before_rr", "before_rr", "", "", "fasta"))
        if dataset_info.mode == "standard":
            contigs.append(("assembly_graph", "assembly_graph", "", "", "fastg"))
            contigs.append(("assembly_graph_with_scaffolds", "assembly_graph_with_scaffolds", "", "", "gfa"))
    return contigs


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('info', metavar='CONFIG_FILE', type=str,  help='teamcity.py info config')
    parser.add_argument("--run_name", "-n", help="output dir custom suffix", type=str)
    parser.add_argument("--spades_path", "-p", help="custom directory to spades.py", type=str)
    parser.add_argument("--no_cfg_and_compilation", dest="cfg_compilation", help="don't copy configs or compile SPAdes even if .info file states otherwise", action='store_false')
    parser.set_defaults(cfg_compilation=True)
    parser.add_argument("--no_contig_archive", dest="contig_archive", help="don't save contigs to common contig archive", action='store_false')
    parser.set_defaults(contig_archive=True)
    parser.add_argument("--contig_name", "-s", help="archive contig name custom suffix", type=str)
    parser.add_argument("--spades_cfg_dir", "-c", help="SPAdes config directory", type=str)
    parser.add_argument("--local_output_dir", "-o", help="use this output dir, override output directory provided in config", type=str)
    args = parser.parse_args()
    return args


# Save meta information about this teamcity.py run
def save_run_info(args, output_dir):
    run_info = open(os.path.join(output_dir, "test_run.info"), "w")
    run_info.write(".info file: " + args.info + "\n");
    if args.run_name:
        run_info.write("run name: " + args.run_name + "\n");
    if args.spades_path:
        run_info.write("path to spades.py: " + str(args.spades_path) + "\n");
    if args.contig_archive:
        run_info.write("save contigs archive: " + str(args.contig_archive) + "\n");
    if args.contig_name:
        run_info.write("contig custom name: " + args.contig_name + "\n");
    if args.spades_cfg_dir:
        run_info.write("spades config direrctory: " + args.spades_cfg_dir + "\n");
    if args.local_output_dir:
        run_info.write("local output dir: " + args.local_output_dir + "\n");
    run_info.close()


# Find path to contigs archive
def get_contigs_storage_dir(args, dataset_info):
    contig_storage_dir = ''
    if 'contig_storage' in dataset_info.__dict__ and args.contig_archive:
        contig_storage_dir = dataset_info.contig_storage
        if args.run_name:
            contig_storage_dir += "_" + args.run_name
    return contig_storage_dir


# Create output folder
def create_output_dir(args, dataset_info):
    #make dirs and remembering history
    output_dir = dataset_info.name
    if 'build_agent' in dataset_info.__dict__:
        output_dir += "_" + dataset_info.build_agent

    #add custom suffix if present
    if args.run_name:
        output_dir += "_" + args.run_name

    #override output dir set by .info config
    if args.local_output_dir:
        output_dir = os.path.join(args.local_output_dir, output_dir)
    else:
        output_dir = os.path.join(dataset_info.output_dir, output_dir)

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    return output_dir


### Pipeline functions ###

# Compile SPAdes
def compile_spades(args, dataset_info, working_dir):
    if not args.cfg_compilation:
        log.log("Forced to use current SPAdes build, will not compile SPAdes");
    elif 'spades_compile' not in dataset_info.__dict__ or dataset_info.spades_compile:
        comp_params = ' '
        if 'compilation_params' in dataset_info.__dict__:
            comp_params = " ".join(dataset_info.compilation_params)

        bin_dir = 'build_spades'
        if not os.path.exists(bin_dir):
            os.makedirs(bin_dir)
        os.chdir(bin_dir)

        #Compilation
        err_code = os.system('cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=' + working_dir + ' ' + os.path.join(working_dir, 'src') + comp_params)
        err_code = err_code | os.system('make -j 16')
        err_code = err_code | os.system('make install')

        os.chdir(working_dir)

        if err_code != 0:
            # Compile from the beginning if failed
            shutil.rmtree('bin', True)
            shutil.rmtree('build_spades', True)
            return os.system('./spades_compile.sh ' + comp_params)
    return 0


#Create SPAdes command line
def make_spades_cmd(args, dataset_info, spades_dir, output_dir):
    #Correct paths in params
    input_file_option_list = ['-1', '-2', '--12', '-s', '--hap' , '--trusted-contigs', '--untrusted-contigs' , '--nanopore' , '--pacbio', '--sanger', '--pe1-1', '--pe1-2', '--pe1-s', '--pe2-1', '--pe2-2', '--pe2-s', '--mp1-1', '--mp1-2', '--mp1-s', '--mp2-1', '--mp2-2', '--mp2-s', '--hqmp1-1', '--hqmp1-2', '--hqmp1-s', '--hqmp2-1', '--hqmp2-2', '--hqmp2-s', '--tslr', '--dataset'];
    spades_params = []
    i = 0
    while i < len(dataset_info.spades_params):
        option = dataset_info.spades_params[i]
        spades_params.append(str(option))
        if i < len(dataset_info.spades_params) - 1 and option in input_file_option_list:
            spades_params.append(os.path.join(dataset_info.dataset_path, str(dataset_info.spades_params[i + 1])))
            i += 1
        i += 1

    if args.spades_cfg_dir:
        spades_params.append("--configs-dir")
        spades_params.append(args.spades_cfg_dir)

    spades_exec = dataset_info.mode + "spades.py"
    if dataset_info.mode in ['standard', 'tru']:
        spades_exec = "spades.py"
    if "--only-assembler" not in spades_params:
        spades_exec += " --disable-gzip-output "

    return os.path.join(spades_dir, spades_exec) + " " + " ".join(spades_params) + " -o " + output_dir


# Check etalon saves using detect_diffs.sh
def check_etalon_saves(dataset_info, working_dir, output_dir):
    exit_code = 0
    log.log("Comparing etalon saves now")
    ecode = os.system(os.path.join(working_dir, "src/test/teamcity/detect_diffs.sh") + " " + output_dir + " " + dataset_info.etalon_saves)
    if ecode != 0:
        log.err("Comparing etalon saves finished abnormally with exit code " + str(ecode))
        exit_code = 12
    return exit_code


# Save contigs to storage
def save_contigs(args, output_dir, contig_storage_dir, contigs, rewrite_latest):
    if contig_storage_dir == '':
        return

    if not os.path.exists(contig_storage_dir):
        os.makedirs(contig_storage_dir)

    name_prefix = datetime.datetime.now().strftime('%Y%m%d-%H%M') + "_"
    log.log("Contigs have prefix " + name_prefix)
    ctg_suffix = ""
    if args.contig_name:
        ctg_suffix = "_" + args.contig_name
        log.log("Contigs have suffix " + ctg_suffix)

    for name, file_name, prefix, opts, ext in contigs:
        saved_ctg_name = name_prefix + name + ctg_suffix + "." + ext
        spades_filename = os.path.join(output_dir, file_name + "." + ext)
        if not os.path.exists(spades_filename):
            log.warn("File " + spades_filename + " do not exist ")
            continue
        shutil.copy(spades_filename, os.path.join(contig_storage_dir, saved_ctg_name))
        log.log(name + " saved to " + os.path.join(contig_storage_dir, saved_ctg_name))

        if rewrite_latest:
            log.log("Creating latest symlink")
            lc =  os.path.join(contig_storage_dir, "latest_" + name + "." + ext)
            if os.path.islink(lc):
                os.remove(lc)
            os.symlink(os.path.join(contig_storage_dir, saved_ctg_name), lc)


# Save QUAST report as artifact
def save_quast_report(contigs, dataset_info, contig_storage_dir, output_dir, artifact_dir, clean_artifacts_dir = True):
    working_dir = os.getcwd()
    if not (dataset_info.mode in ("standard", "tru", "meta", "plasmid")):
        return

    if clean_artifacts_dir:
        shutil.rmtree(artifact_dir, True)

    if not os.path.exists(artifact_dir):
        os.makedirs(artifact_dir)

    log.log("======= Saving report of to " + artifact_dir + " =======")
    for name, file_name, prefix, opts, ext in contigs:
        if ext != "fasta":
            continue
        quast_output_dir = os.path.join(output_dir, "QUAST_RESULTS_CMP_" + name.upper())

        if not os.path.exists(quast_output_dir):
            if not ('quast_params' in dataset_info.__dict__ and dataset_info.quast_params):
                log.err('QUAST params are not set')
                return

            fn = os.path.join(output_dir, file_name + "." + ext)
            if not os.path.exists(fn):
                log.err('File for analysis is no found: ' + fn)
                continue

            latest_ctg = os.path.join(contig_storage_dir, "latest_" + name + "." + ext)
            flist = [fn]
            if not os.path.exists(fn):
                flist.append(latest_ctg)

            run_quast(dataset_info, flist, quast_output_dir, opts)

        compressed_report = os.path.join(artifact_dir, name + "_report.zip")
        if os.path.exists(compressed_report):
            os.remove(compressed_report)

        log.log("Saving report of " + name + " to " + compressed_report)
        os.chdir(quast_output_dir)
        os.system("zip -9r " + compressed_report + " * > /dev/null")
        os.chdir(working_dir)

