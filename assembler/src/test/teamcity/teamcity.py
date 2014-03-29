#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
        log_str += metric + " = " + str(result_map[metric])

        if metric in limit_map:
            if limit_map[metric][1]:
                if result_map[metric] < limit_map[metric][0]:
                    print(metric + " = " + str(result_map[metric]) + " is less than expected: " + str(limit_map[metric][0]))
                    log_str += " (bad)"
                    res = -1
                else:
                    print(metric + " = " + str(result_map[metric]) + " >= " + str(limit_map[metric][0]) + " (OK)")

            else:
                if result_map[metric] > limit_map[metric][0]:
                    print(metric + " = " + str(result_map[metric]) + " is higher than expected: " + str(limit_map[metric][0]))
                    log_str += " (bad)"
                    res = -1
                else:
                    print(metric + " = " + str(result_map[metric]) + " <= " + str(limit_map[metric][0]) + " (OK)")
        else:
            print(metric + " = " + str(result_map[metric]))


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
    sys.exit(1)

dataset_path, dataset_info = load_info(sys.argv[1])

#clean
if ('prepare_cfg' not in dataset_info.__dict__ or dataset_info.prepare_cfg) and ('spades_compile' not in dataset_info.__dict__ or dataset_info.spades_compile):
    shutil.rmtree('bin', True)
    shutil.rmtree('build', True)
    shutil.rmtree('build_spades', True)

#prepare cfg
if 'prepare_cfg' not in dataset_info.__dict__ or dataset_info.prepare_cfg:
    ecode = os.system('./prepare_cfg')
    if ecode != 0:
        print("Preparing configuration files finished abnormally with exit code " + str(ecode))
        sys.exit(2)

#compile
if 'spades_compile' not in dataset_info.__dict__ or dataset_info.spades_compile:
    comp_params = ' '
    if 'compilation_params' in dataset_info.__dict__:
        comp_params = " ".join(dataset_info.compilation_params)

    ecode = os.system('./spades_compile.sh ' + comp_params)
    if ecode != 0:
        print("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)

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
    if i < len(dataset_info.spades_params) - 1 and (option == '-1' or option == '-2' or option == '--12' or option == '-s' or option == '--hap'  or option == '--trusted-contigs' or option == '--untrusted-contigs'  or option == '--pacbio' or option == '--sanger'):
        spades_params.append(os.path.join(dataset_path, str(dataset_info.spades_params[i + 1])))
        i += 1
    i += 1

spades_cmd = ""
if 'dipspades' in dataset_info.__dict__ and dataset_info.dipspades:
    spades_cmd = "./src/spades_pipeline/dipspades_logic.py " + " ".join(spades_params) + " -o " + output_dir
else:
    spades_cmd = "./spades.py --disable-gzip-output " + " ".join(spades_params) + " -o " + output_dir

#run spades
ecode = os.system(spades_cmd)
if ecode != 0:
    print("SPAdes finished abnormally with exit code " + str(ecode))
    write_log(history_log, "", output_dir, dataset_info)
    os.system("chmod -R 777 " + output_dir)
    sys.exit(4)

new_log = ''
exit_code = 0

#reads quality
if 'reads_quality_params' in dataset_info.__dict__:
    corrected_reads_dataset = os.path.join(output_dir, "corrected/corrected.yaml")

    if not os.path.exists(corrected_reads_dataset):
        print("Corrected reads were not detected in " + corrected_reads_dataset)
        exit_code = 5
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
            sys.exit(6)

        limit_map = {}
        if 'assess' in dataset_info.__dict__ and dataset_info.assess:
            if 'min_genome_mapped' in dataset_info.__dict__:
                limit_map["Genome mapped"] = (float(dataset_info.min_genome_mapped), True)
            if 'min_aligned' in dataset_info.__dict__:
                limit_map["Uniquely aligned reads"] = (float(dataset_info.min_aligned), True)
            
        result = assess_reads(os.path.join(rq_output_dir, "report.horizontal.tsv"), limit_map)
        if result[0] != 0:
            exit_code = 7
        new_log += result[1]


#QUAST
quast_cmd = ""
if 'quast_params' in dataset_info.__dict__:
    contigs = os.path.join(output_dir, "contigs.fasta")
    if 'dipspades' in dataset_info.__dict__ and dataset_info.dipspades:
        contigs = os.path.join(output_dir, "consensus_contigs.fasta")

    if not os.path.exists(contigs):
        print("No contigs were found in " + output_dir)
        exit_code = 8
    else:
        quast_params = []
        if dataset_info.quast_params:
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
            sys.exit(9)

        limit_map = {}
        if 'assess' in dataset_info.__dict__ and dataset_info.assess:
            print("Assessing QUAST results...")
            limit_map = {}
            if 'min_n50' in dataset_info.__dict__:
                limit_map["N50"] = (int(dataset_info.min_n50), True)
            if 'max_mis' in dataset_info.__dict__:
                limit_map["Misassemblies"] = (int(dataset_info.max_mis), False)
            if 'min_genome_mapped' in dataset_info.__dict__:
                limit_map["Genome mapped"] = (float(dataset_info.min_genome_mapped), True)
            if 'min_genes' in dataset_info.__dict__:
                limit_map["Genes"] = (int(dataset_info.min_genes), True)
            if 'max_indels' in dataset_info.__dict__:
                limit_map["Indels"] = (float(dataset_info.max_indels), False)
            if 'max_subs' in dataset_info.__dict__:
                limit_map["Mismatches"] = (float(dataset_info.max_subs), False)

        result = assess_quast(os.path.join(quast_output_dir, "transposed_report.tsv"), limit_map, "contigs")
        if result[0] != 0:
            exit_code = 10
        new_log += result[1]

        #SCAFFOLDS
        scafs = os.path.join(output_dir, "scaffolds.fasta")
        if os.path.exists(scafs):
            quast_output_scaf_dir = os.path.join(output_dir, "QUAST_RESULTS_SCAF")
            if os.system(quast_cmd + " -o " + quast_output_scaf_dir + " " + scafs) != 0:
                print("Failed to estimate scaffolds")
                if 'sc_assess' in dataset_info.__dict__ and dataset_info.sc_assess:
                    exit_code = 11
            else:
                sc_limit_map = {}
                if 'sc_assess' in dataset_info.__dict__ and dataset_info.sc_assess:
                    print("Assessing QUAST results for scaffolds...")
                    if 'sc_min_n50' in dataset_info.__dict__:
                        sc_limit_map["N50"] = (int(dataset_info.sc_min_n50), True)
                    if 'sc_max_mis' in dataset_info.__dict__:
                        sc_limit_map["Misassemblies"] = (int(dataset_info.sc_max_mis), False)
                    if 'sc_min_genome_mapped' in dataset_info.__dict__:
                        sc_limit_map["Genome mapped"] = (float(dataset_info.sc_min_genome_mapped), True)
                    if 'sc_min_genes' in dataset_info.__dict__:
                        sc_limit_map["Genes"] = (int(dataset_info.sc_min_genes), True)
                    if 'sc_max_indels' in dataset_info.__dict__:
                        sc_limit_map["Indels"] = (float(dataset_info.sc_max_indels), False)
                    if 'sc_max_subs' in dataset_info.__dict__:
                        sc_limit_map["Mismatches"] = (float(dataset_info.sc_max_subs), False)
                result = assess_quast(os.path.join(quast_output_scaf_dir, "transposed_report.tsv"), sc_limit_map, "scaffolds")
                if result[0] != 0:
                    exit_code = 11
                new_log += result[1]


#etalon saves
if 'etalon_saves' in dataset_info.__dict__:
    print("Comparing etalon saves now")
    ecode = os.system("./src/test/teamcity/detect_diffs.sh " + output_dir + " " + dataset_info.etalon_saves)
    if ecode != 0:
        print("Comparing etalon saves finished abnormally with exit code " + str(ecode))
        exit_code = 12


#writing log
write_log(history_log, new_log, output_dir, dataset_info)


#saving contigs
if 'contig_storage' in dataset_info.__dict__:
    contig_dir = dataset_info.contig_storage
    quast_contig_dir = os.path.join(contig_dir, "quast_results")
    if not os.path.exists(quast_contig_dir):
        os.makedirs(quast_contig_dir)

    name_prefix = datetime.datetime.now().strftime('%Y%m%d-%H%M')
    if len(sys.argv) == 3:
        name_prefix += "_" + sys.argv[2]
    print("Contigs have prefix " + name_prefix)

    if 'dipspades' in dataset_info.__dict__ and dataset_info.dipspades:
        shutil.copy(os.path.join(output_dir, "consensus_contigs.fasta"), os.path.join(contig_dir, name_prefix + ".fasta"))
    else:
        shutil.copy(os.path.join(output_dir, "contigs.fasta"), os.path.join(contig_dir, name_prefix + ".fasta"))
    print("Contigs saved to " + os.path.join(contig_dir, name_prefix + ".fasta"))

    scafs = os.path.join(output_dir, "scaffolds.fasta")
    if os.path.exists(scafs):
        shutil.copy(scafs, os.path.join(contig_dir, name_prefix + "_scafs.fasta"))
        print("Scaffolds saved to " + os.path.join(contig_dir, name_prefix + "_scafs.fasta"))

    before_rr = os.path.join(output_dir, "before_rr.fasta")
    if os.path.exists(before_rr):
        shutil.copy(before_rr, os.path.join(contig_dir, name_prefix + "_before_rr.fasta"))
        print("Contigs before resolve saved to " + os.path.join(contig_dir, name_prefix + "_before_rr.fasta"))

    before_corr = os.path.join(output_dir, "assembled_contigs.fasta")
    if os.path.exists(before_corr):
        shutil.copy(before_corr, os.path.join(contig_dir, name_prefix + "_before_corr.fasta"))

    scaff_before_corr = os.path.join(output_dir, "assembled_scaffolds.fasta")
    if os.path.exists(scaff_before_corr):
        shutil.copy(scaff_before_corr, os.path.join(contig_dir, name_prefix + "_scaff_before_corr.fasta"))

    boken_scaff = os.path.join(output_dir, "broken_scaffolds.fasta")
    if os.path.exists(boken_scaff):
        shutil.copy(boken_scaff, os.path.join(contig_dir, name_prefix + "_scaff_broken.fasta"))


#    if quast_cmd != "":
#        import glob
#        sys.path.append('./src/tools/contig_analysis/')
#        import compare_fasta
#        print("Running quast for all saved contigs now...")
#        os.system(quast_cmd + " -o " + quast_contig_dir + " " + os.path.join(contig_dir, "*.fasta") + " > " + os.path.join(contig_dir, "quast.log") + " 2> " + os.path.join(contig_dir, "quast.err"))
#        print("Done")

os.system("chmod -R 777 " + output_dir)
sys.exit(exit_code)

