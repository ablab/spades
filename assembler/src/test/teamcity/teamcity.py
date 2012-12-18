#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#script for testing SPAdes
#provide a directory with dataset.info file, or complete path to .info file


import sys
import os
import shutil
import getopt
import re
import datetime

sys.path.append('./src/spades_pipeline/')
import process_cfg


def load_info(dataset_path):
    if not os.path.exists(dataset_path):
        print(dataset_path + " can not be found")
        sys.exit(-2)

    if os.path.isdir(dataset_path):
        dataset_path = os.path.join(dataset_path, "dataset.info")

    return os.path.split(dataset_path)[0], process_cfg.load_config_from_file(dataset_path)


def read_log(output_dir, dataset_info):
    history_log = ""
    if 'history_log' in dataset_info.__dict__:
        log_file_name = dataset_info.history_log
        if log_file_name[0] != '/':
            log_file_name = os.path.join(output_dir, log_file_name)

        if os.path.exists(log_file_name) and os.path.isfile(log_file_name):
            history_log = ""
            log_file = open(log_file_name, "r")
            for line in log_file:
                history_log += line
            log_file.close()

    return history_log


def write_log(history_log, new_log, output_dir, dataset_info):
    if 'history_log' in dataset_info.__dict__:
        log_file_name = dataset_info.history_log
        log_file_name = os.path.join(output_dir, log_file_name)

        log_file = open(log_file_name, "w")
        log_file.write(history_log)
        log_file.write(new_log)
        log_file.close()


def assess_map(result_map, limit_map):
    res = 0
    log_str = ""

    for metric in sorted(result_map.keys()):
        print(metric + " = " + str(result_map[metric]))
        log_str += metric + " = " + str(result_map[metric])

        if metric in limit_map:
            if limit_map[metric][1] and result_map[metric] < limit_map[metric][0]:
                print(metric + " is less than expected: " + str(limit_map[metric][0]))
                log_str += " (bad)"
                res -= 1
            elif not limit_map[metric][1] and result_map[metric] > limit_map[metric][0]:
                print(metric + " is higher than expected: " + str(limit_map[metric][0]))
                log_str += " (bad)"
                res -= 1

        log_str += "; "

    return res, log_str
                

def assess_reads(report, limit_map):
    f = open(report, 'r')
    columns = map(lambda s: s.strip(), f.readline().split('\t'))
    values = map(lambda s: s.strip(), f.readline().split('\t'))
    f.close()

    result_map = {}
    if "Aligned read" in columns:
        aligned = values[columns.index("Aligned reads")]
        aligned = re.search('\((.+)%\)', aligned).group(1)
        result_map["Aligned reads"] = float(aligned)
    if "Genome mapped (%)" in columns:
        result_map["Genome mapped"] = float(values[columns.index("Genome mapped (%)")])

    log_str = "Read quality " + str(datetime.datetime.now().strftime('%d.%m.%Y %H:%M')) + "\n"

    result = assess_map(result_map, limit_map)

    log_str += result[1] + ".\n"
    return result[0], log_str


def assess_quast(report, limit_map, name = ""):
    f = open(report, 'r')
    columns = map(lambda s: s.strip(), f.readline().split('\t'))
    values = map(lambda s: s.strip(), f.readline().split('\t'))
    f.close()

    result_map = {}
    if "N50" in columns:
        result_map["N50"] = int(values[columns.index("N50")])
    if "# misassemblies" in columns:
        result_map["Misassemblies"] = int(values[columns.index("# misassemblies")])
    if "# genes" in columns:
        result_map["Genes"] = int(values[columns.index("# genes")].split('+')[0])
    if "Genome fraction (%)" in columns:
        result_map["Genome mapped"] = float(values[columns.index("Genome fraction (%)")])
    if "# mismatches per 100 kbp" in columns:
        result_map["Mismatches"] = float(values[columns.index("# mismatches per 100 kbp")])
    if "# indels per 100 kbp" in columns:
        result_map["Indels"] = float(values[columns.index("# indels per 100 kbp")])

    log_str = "Quast " + str(datetime.datetime.now().strftime('%d.%m.%Y %H:%M')) + " " + name + "\n"

    print("Assessing " + name)
    result = assess_map(result_map, limit_map)

    log_str += result[1] + ".\n"
    return result[0], log_str


def filter_file_list(file_list):
    current = 0
    result_list = [file_list[current]]    
    current_set = compare_fasta.read_fasta(file_list[current])

    i = 1
    while i < len(file_list):
        if not compare_contigs(current_set, read_fasta(file_list[i])):
            current = i
            result_list.append(file_list[current])
            current_set = compare_fasta.read_fasta(file_list[current])
        i += 1
    
    return result_list


### main ###

if len(sys.argv) < 2:
    print("Data set info is not provided")
    sys.exit(-1)

dataset_path, dataset_info = load_info(sys.argv[1])

#prepare cfg
if 'prepare_cfg' not in dataset_info.__dict__ or dataset_info.prepare_cfg:
    ecode = os.system('./prepare_cfg')
    if ecode != 0:
        print("Preparing configuration files finished abnormally with exit code " + str(ecode))
        sys.exit(ecode)


#compile
if 'spades_compile' not in dataset_info.__dict__ or dataset_info.spades_compile:
    ecode = os.system('./spades_compile.sh')
    if ecode != 0:
        print("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(ecode)

#make dirs and remembering history
spades_output_dir_name = dataset_info.name
if 'build_agent' in dataset_info.__dict__:
    spades_output_dir_name += "_" + dataset_info.build_agent
output_dir = os.path.join(dataset_info.output_dir, spades_output_dir_name)

history_log = read_log(output_dir, dataset_info)
if history_log != "":
    print("Quality log found, going to append")
    
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.makedirs(output_dir)
os.system("chmod -R 777 " + output_dir)

#make correct files to files
spades_params = []
i = 0
while i < len(dataset_info.spades_params):
    option = dataset_info.spades_params[i]
    spades_params.append(str(option))
    if i < len(dataset_info.spades_params) - 1 and (option == '-1' or option == '-2' or option == '--12' or option == '-s'):
        spades_params.append(os.path.join(dataset_path, str(dataset_info.spades_params[i + 1])))
        i += 1
    i += 1

spades_cmd = "./spades.py --disable-gzip-output " + " ".join(spades_params) + " -o " + output_dir

#run spades
ecode = os.system(spades_cmd)
if ecode != 0:
    print("SPAdes finished abnormally with exit code " + str(ecode))
    write_log(history_log, "", output_dir, dataset_info)
    os.system("chmod -R 777 " + output_dir)
    sys.exit(ecode)

new_log = ''
exit_code = 0

#reads quality
if 'reads_quality_params' in dataset_info.__dict__:
    corrected_reads_dataset = os.path.join(output_dir, "corrected/dataset.info")

    if not os.path.exists(corrected_reads_dataset):
        print("Corrected reads were not detected in " + corrected_reads_dataset)
        exit_code = -3
    else:
        rq_params = []
        i = 0
        while i < len(dataset_info.reads_quality_params):
            option = dataset_info.reads_quality_params[i]
            rq_params.append(str(option))
            if i < len(dataset_info.reads_quality_params) - 1 and option == '-r':
                rq_params.append(os.path.join(dataset_path, str(dataset_info.reads_quality_params[i + 1])))
                i += 1
            i += 1

        rq_output_dir = os.path.join(output_dir, "RQ_RESULTS")
        rq_cmd = "./src/tools/reads_utils/reads_quality.py " + " ".join(rq_params) + " -o " + rq_output_dir + " " + corrected_reads_dataset
        ecode = os.system(rq_cmd)
        if ecode != 0:
            print("Reads quality tool finished abnormally with exit code " + str(ecode))
            write_log(history_log, "", output_dir, dataset_info)
            os.system("chmod -R 777 " + output_dir)
            sys.exit(ecode)

        limit_map = {}
        if 'assess' in dataset_info.__dict__ and dataset_info.assess:
            if 'min_genome_mapped' in dataset_info.__dict__:
                limit_map["Genome mapped"] = (dataset_info.min_genome_mapped, True)
            if 'min_aligned' in dataset_info.__dict__:
                limit_map["Aligned reads"] = (dataset_info.min_aligned, True)
            
        result = assess_reads(os.path.join(rq_output_dir, "report.horizontal.tsv"), limit_map)
        exit_code += result[0]
        new_log += result[1]


#QUAST
quast_cmd = ""
if 'quast_params' in dataset_info.__dict__:
    contigs = os.path.join(output_dir, "contigs.fasta")

    if not os.path.exists(contigs):
        print("No contigs were found in " + output_dir)
        exit_code = -4
    else:
        quast_params = []
        i = 0
        while i < len(dataset_info.quast_params):
            option = dataset_info.quast_params[i]
            quast_params.append(str(option))
            if i < len(dataset_info.quast_params) - 1 and (option == '-R' or option == '-G' or option == '-O'):
                quast_params.append(os.path.join(dataset_path, str(dataset_info.quast_params[i + 1])))
                i += 1
            i += 1

        #CONTIGS
        quast_output_dir = os.path.join(output_dir, "QUAST_RESULTS")
        quast_cmd = os.path.join(dataset_info.quast_dir, "quast.py") + " " + " ".join(quast_params)
        ecode = os.system(quast_cmd + " -o " + quast_output_dir + " " + contigs)
        if ecode != 0:
            print("QUAST finished abnormally with exit code " + str(ecode))
            write_log(history_log, "", output_dir, dataset_info)
            os.system("chmod -R 777 " + output_dir)
            sys.exit(ecode)

        limit_map = {}
        if 'assess' in dataset_info.__dict__ and dataset_info.assess:
            limit_map = {}
            if 'min_n50' in dataset_info.__dict__:
                limit_map["N50"] = (dataset_info.min_n50, True)
            if 'max misassemblies' in dataset_info.__dict__:
                limit_map["Misassemblies"] = (dataset_info.max_misassemblies, False)
            if 'min_genome_mapped' in dataset_info.__dict__:
                limit_map["Genome mapped"] = (dataset_info.min_genome_mapped, True)
            if 'min genes ' in dataset_info.__dict__:
                limit_map["Genes"] = (dataset_info.min_genes, True)
            if 'max indels' in dataset_info.__dict__:
                limit_map["Mismatches"] = (dataset_info.max_indels, False)
            if 'max subs' in dataset_info.__dict__:
                limit_map["Mismatches"] = (dataset_info.max_subs, False)
            
        result = assess_quast(os.path.join(quast_output_dir, "transposed_report.tsv"), limit_map, "contigs")
        exit_code += result[0]
        new_log += result[1]

        #SCAFFOLDS
        scafs = os.path.join(output_dir, "scaffolds.fasta")
        if os.path.exists(scafs):
            quast_output_scaf_dir = os.path.join(output_dir, "QUAST_RESULTS_SCAF")
            quast_cmd = os.path.join(dataset_info.quast_dir, "quast.py") + " " + " ".join(quast_params) + " -o " + quast_output_scaf_dir + " " + scafs
            if os.system(quast_cmd) != 0:
                print("Failed to estimate scaffolds")
            else:
                result = assess_quast(os.path.join(quast_output_scaf_dir, "transposed_report.tsv"), {}, "scaffolds")
                new_log += result[1]


#etalon saves
if 'etalon_saves' in dataset_info.__dict__:
    ecode = os.system("./src/test/teamcity/detect_diffs.sh " + output_dir + " " + dataset_info.etalon_saves)
    if ecode != 0:
        print("Comparing etalon saves finished abnormally with exit code " + str(ecode))
        exit_code = ecode


#writing log
write_log(history_log, new_log, output_dir, dataset_info)


#saving contigs
if 'contig_storage' in dataset_info.__dict__:
    contig_dir = dataset_info.contig_storage
    quast_contig_dir = os.path.join(contig_dir, "quast_results")
    if not os.path.exists(quast_contig_dir):
        os.makedirs(quast_contig_dir)

    name_prefix = datetime.datetime.now().strftime('%Y%m%d-%H%M')
    shutil.copy(os.path.join(output_dir, "contigs.fasta"), os.path.join(contig_dir, name_prefix + ".fasta"))

    scafs = os.path.join(output_dir, "scaffolds.fasta")
    if os.path.exists(scafs):
        shutil.copy(scafs, os.path.join(contig_dir, name_prefix + "_scafs.fasta"))

    if quast_cmd != "":
        import glob
#        sys.path.append('./src/tools/contig_analysis/')
#        import compare_fasta
        print("Running quast for all saved contigs now...")
        os.system(quast_cmd + " -o " + quast_contig_dir + " " + os.path.join(contig_dir, "*.fasta") + " > " + os.path.join(contig_dir, "quast.log") + " 2> " + os.path.join(contig_dir, "quast.err"))
        print("Done")

os.system("chmod -R 777 " + output_dir)
sys.exit(exit_code)

