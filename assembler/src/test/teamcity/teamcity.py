#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
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
from traceback import print_exc

sys.path.append('./src/spades_pipeline/')
import process_cfg


class Log:

    text = ""


    def log(self, s):
        self.text += s + "\n"    


    def print_log(self):
        print(self.text)


log = Log()


def load_info(dataset_path):
    if not os.path.exists(dataset_path):
        print(dataset_path + " can not be found")
        sys.exit(-2)

    if os.path.isdir(dataset_path):
        dataset_path = os.path.join(dataset_path, "dataset.info")

    info = process_cfg.load_config_from_file(dataset_path)
    info.__dict__["dataset_path"] = os.path.split(dataset_path)[0]
    return info

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

        if metric in limit_map and len(limit_map[metric]) > 0:
            for constraint in limit_map[metric]:
                if constraint[1]:
                    if result_map[metric] < constraint[0]:
                        log.log(metric + " = " + str(result_map[metric]) + " is less than expected: " + str(constraint[0]))
                        log_str += " (bad)"
                        res = -1
                    else:
                        log.log(metric + " = " + str(result_map[metric]) + " >= " + str(constraint[0]) + " (OK)")

                else:
                    if result_map[metric] > constraint[0]:
                        log.log(metric + " = " + str(result_map[metric]) + " is higher than expected: " + str(constraint[0]))
                        log_str += " (bad)"
                        res = -1
                    else:
                        log.log(metric + " = " + str(result_map[metric]) + " <= " + str(constraint[0]) + " (OK)")
        else:
            log.log(metric + " = " + str(result_map[metric]))


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
        result_map["N50"] = float(values[columns.index("N50")])
    if "# misassemblies" in columns:
        result_map["Misassemblies"] = float(values[columns.index("# misassemblies")])
    if "# genes" in columns:
        result_map["Genes"] = float(values[columns.index("# genes")].split('+')[0])
    if "Genome fraction (%)" in columns:
        result_map["Genome mapped"] = float(values[columns.index("Genome fraction (%)")])
    if "# mismatches per 100 kbp" in columns:
        result_map["Mismatches"] = float(values[columns.index("# mismatches per 100 kbp")])
    if "# indels per 100 kbp" in columns:
        result_map["Indels"] = float(values[columns.index("# indels per 100 kbp")])

    log_str = "Quast " + str(datetime.datetime.now().strftime('%d.%m.%Y %H:%M')) + " " + name + "\n"

    log.log("Assessing " + name)
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


def run_quast(dataset_info, contigs, quast_output_dir):
    if not reduce(lambda x, y: os.path.exists(y) and x, contigs, True):
        print("No contigs were found in " + output_dir)
        return 8
    else:
        cmd = "quast.py"
        if 'meta' in dataset_info.__dict__ and dataset_info.meta:
            cmd = "metaquast.py"
        quast_params = []
        if dataset_info.quast_params:
            i = 0
            while i < len(dataset_info.quast_params):
                option = dataset_info.quast_params[i]
                quast_params.append(str(option))
                if i < len(dataset_info.quast_params) - 1 and (option == '-R' or option == '-G' or option == '-O'):
                    quast_params.append(os.path.join(dataset_info.dataset_path, str(dataset_info.quast_params[i + 1])))
                    i += 1
                i += 1

        quast_cmd = os.path.join(dataset_info.quast_dir, cmd) + " " + " ".join(quast_params)
        ecode = os.system(quast_cmd + " -o " + quast_output_dir + " " + " ".join(contigs) )
        if ecode != 0:
            print("QUAST finished abnormally with exit code " + str(ecode))
            return 9
    return 0


def construct_map(dataset_info, prefix):
    limit_map = {}
    params = map(lambda x: (prefix + x[0],x[1],x[2]),
                                       [('min_n50', "N50", True) , 
                                        ('max_n50', "N50", False),
                                        ('max_mis', "Misassemblies", False),
                                        ('min_mis', "Misassemblies", True),
                                        ('min_genome_mapped', "Genome mapped", True),
                                        ('min_genes', "Genes", True),       
                                        ('max_indels', "Indels", False), 
                                        ('max_subs', "Mismatches", False)])

    if prefix + 'assess' in dataset_info.__dict__ and dataset_info.__dict__[prefix + 'assess']:
        print("Assessing QUAST results...")
        limit_map = {"N50":[], "Misassemblies":[], "Genome mapped":[], "Genes":[], "Indels":[], "Mismatches":[]}
        
        for p in params:
            if p[0] in dataset_info.__dict__:
                limit_map[p[1]].append((float(dataset_info.__dict__[p[0]]), p[2]))
    return limit_map


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

    # main part of plantagora output
    prev_line = ""
    coords = {}
    currend_id = ""
    while not line.startswith("Analyzing coverage..."):
        line = infile.readline()
        if line.startswith("CONTIG:"):
            currend_id = line.split()[1].strip()
        elif (line.find("Extensive misassembly") != -1):
            line2 = infile.readline()
            mis = parse_misassembly(prev_line, line2)
            if not mis:
                print("Failed to parse misassembly")
            elif mis not in coords:
                coords[mis] = [currend_id]
            else:
                coords[mis].append(currend_id)
            prev_line = line2
            continue
        prev_line = line
               
    infile.close()
    return coords

def quast_run_and_assess(dataset_info, fn, output_dir, name, prefix, special_exit_code):
    if os.path.exists(fn):
        print("Processing " + fn)
        qcode = run_quast(dataset_info, [fn], output_dir)
        if qcode != 0:
            print("Failed to estimate!")
            if (prefix + 'assess') in dataset_info.__dict__ and dataset_info.__dict__[prefix + 'assess']:
                return special_exit_code, new_log
            return qcode, ""
     
        limit_map = construct_map(dataset_info, prefix)
        report_path = output_dir
        if 'meta' in dataset_info.__dict__ and dataset_info.meta:
            report_path = os.path.join(report_path, "combined_quast_output")
        report_path = os.path.join(report_path, "transposed_report.tsv")

        result = assess_quast(report_path, limit_map, name)

        if result[0] != 0:
            return special_exit_code, result[1]

        return 0, result[1]
    else:
        print("File not found " + fn)
        return 8, ""

def quast_analysis(dataset_info, folder):
    exit_code = 0
    new_log = ''
    contigs = "contigs"
    if 'dipspades' in dataset_info.__dict__ and dataset_info.dipspades:
        contigs = "consensus_contigs"
    elif 'truseq' in dataset_info.__dict__ and dataset_info.truseq:
        contigs = "truseq_long_reads"

    log.log("======= CONTIG SUMMARY =======")
    qcode, qlog = quast_run_and_assess(dataset_info, os.path.join(folder, contigs + ".fasta"), 
                    os.path.join(folder, "QUAST_RESULTS"), contigs, "", 10)
    new_log += qlog
    if qcode != 0:
        print("Contig analysis exit with code ", qcode)
        exit_code = qcode

    log.log("======= SCAFFOLD SUMMARY =======")
    qcode, qlog = quast_run_and_assess(dataset_info, os.path.join(folder, "scaffolds.fasta"), 
                    os.path.join(folder, "QUAST_RESULTS_SCAF"), "scaffolds", "sc_", 11)
    new_log += qlog
    if qcode != 0:
        print("Scaffold analysis exit with code ", qcode)
        exit_code = qcode

    if os.path.exists(os.path.join(folder, "first_pe_contigs.fasta")):
        log.log("======= PRELIMINARY SUMMARY =======")
        qcode, qlog = quast_run_and_assess(dataset_info, os.path.join(folder, "first_pe_contigs.fasta"), 
                        os.path.join(folder, "QUAST_RESULTS_PRELIM"), "preliminary", "prelim_", 25)
        new_log += qlog
        if qcode != 0:
            print("Preliminary scaffold analysis exit with code ", qcode)
            exit_code = qcode

    return exit_code, new_log

def cmp_misassemblies(quast_output_dir, old_ctgs, new_ctgs):
    log.log("Comparing misassemblies in " + old_ctgs + " and " + new_ctgs)
    old_pos = find_mis_positions(os.path.join(quast_output_dir, "contigs_reports/contigs_report_" + old_ctgs + ".stdout"))
    new_pos = find_mis_positions(os.path.join(quast_output_dir, "contigs_reports/contigs_report_" + new_ctgs + ".stdout"))

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
        

### main ###
try:
    if len(sys.argv) < 2:
        print("Data set info is not provided")
        sys.exit(1)

    new_log = ''
    exit_code = 0

    dataset_info = load_info(sys.argv[1])
    working_dir = os.getcwd()

    #make dirs and remembering history
    output_dir = dataset_info.name
    if 'build_agent' in dataset_info.__dict__:
        output_dir += "_" + dataset_info.build_agent

    if len(sys.argv) > 2:
        output_dir += "_" + sys.argv[2]

    output_dir = os.path.join(dataset_info.output_dir, output_dir)

    history_log = read_log(output_dir, dataset_info)
    if history_log != "":
        print("Quality log found, going to append")
        
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    os.system("chmod -R 777 " + output_dir)

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

    #make correct files to files
    spades_params = []
    i = 0
    while i < len(dataset_info.spades_params):
        option = dataset_info.spades_params[i]
        spades_params.append(str(option))
        if i < len(dataset_info.spades_params) - 1 and (option == '-1' or option == '-2' or option == '--12' or option == '-s' or option == '--hap'  or option == '--trusted-contigs' or option == '--untrusted-contigs'  or option == '--pacbio' or option == '--sanger' or option == '--pe1-1' or option == '--pe1-2' or option == '--pe1-s' or option == '--pe2-1' or option == '--pe2-2' or option == '--pe2-s' or option == '--mp1-1' or option == '--mp1-2' or option == '--mp1-s' or option == '--mp2-1' or option == '--mp2-2' or option == '--mp2-s' or option == '--hqmp1-1' or option == '--hqmp1-2' or option == '--hqmp1-s' or option == '--hqmp2-1' or option == '--hqmp2-2' or option == '--hqmp2-s' or option == '--dataset' ):
            spades_params.append(os.path.join(dataset_info.dataset_path, str(dataset_info.spades_params[i + 1])))
            i += 1
        i += 1

    if len(sys.argv) > 3:
        spades_params.append("--configs-dir")
        spades_params.append(sys.argv[3])

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
                    rq_params.append(os.path.join(dataset_info.dataset_path, str(dataset_info.reads_quality_params[i + 1])))
                    i += 1
                i += 1

            rq_output_dir = os.path.join(output_dir, "RQ_RESULTS")
            rq_cmd = "./src/tools/reads_utils/reads_quality.py " + " ".join(rq_params) + " -o " + rq_output_dir + " " + corrected_reads_dataset
            ecode = os.system(rq_cmd)
            if ecode != 0:
                print("Reads quality tool finished abnormally with exit code " + str(ecode))
                write_log(history_log, "", output_dir, dataset_info)
                os.system("chmod -R 777 " + output_dir)
                exit_code = 6
            else:
                limit_map = {}
                if 'assess' in dataset_info.__dict__ and dataset_info.assess:
                    if 'min_genome_mapped' in dataset_info.__dict__:
                        limit_map["Genome mapped"] = [(float(dataset_info.min_genome_mapped), True)]
                    if 'min_aligned' in dataset_info.__dict__:
                        limit_map["Uniquely aligned reads"] = [(float(dataset_info.min_aligned), True)]
                    
                result = assess_reads(os.path.join(rq_output_dir, "report.horizontal.tsv"), limit_map)
                if result[0] != 0:
                    exit_code = 7
                new_log += result[1]


    #QUAST
    if 'quast_params' in dataset_info.__dict__:
        ecode, qlog = quast_analysis(dataset_info, output_dir)
        new_log += qlog
        if ecode != 0:
            print("Quast analysis finished abnormally with exit code " + str(ecode))
            exit_code = ecode
                
    #etalon saves
    if 'etalon_saves' in dataset_info.__dict__:
        print("Comparing etalon saves now")
        ecode = os.system("./src/test/teamcity/detect_diffs.sh " + output_dir + " " + dataset_info.etalon_saves)
        if ecode != 0:
            print("Comparing etalon saves finished abnormally with exit code " + str(ecode))
            exit_code = 12

    contig_dir = ''
    if 'contig_storage' in dataset_info.__dict__:
        contig_dir = dataset_info.contig_storage

        if len(sys.argv) > 2:
            contig_dir += "_" + sys.argv[2]

    #comparing misassemblies
    rewrite_latest = True
    latest_found = True
    if contig_dir != '' and 'quast_params' in dataset_info.__dict__ and '-R' in dataset_info.quast_params and 'assess' in dataset_info.__dict__ and dataset_info.assess:

        latest_ctg = os.path.join(contig_dir, "latest_contigs.fasta")

        if not os.path.exists(latest_ctg):
            import glob
            prev_contigs = sorted(glob.glob(os.path.join(contig_dir, "*_ctg.fasta")))

            if len(prev_contigs) == 0:
                print("No previous contigs detected")
                latest_found = False
            else:
                os.chdir(contig_dir)
                os.symlink(os.path.basename(prev_contigs[-1]), "latest_contigs.fasta")
                os.chdir(working_dir)
        
        if latest_found:
            contigs = "contigs.fasta"
            if 'dipspades' in dataset_info.__dict__ and dataset_info.dipspades:
                contigs = "consensus_contigs.fasta"

            quast_output_dir = os.path.join(output_dir, "QUAST_RESULTS_CMP")
            qcode = run_quast(dataset_info, [os.path.join(output_dir, contigs), latest_ctg], quast_output_dir)

            #compare_misassemblies
            log.log("======= CONTIG COMPARISON =======")
            if not cmp_misassemblies(quast_output_dir, "latest_contigs", os.path.splitext(contigs)[0]):
                rewrite_latest = False
                #exit_code = 13


    if contig_dir != '' and 'quast_params' in dataset_info.__dict__ and '-R' in dataset_info.quast_params and 'sc_assess' in dataset_info.__dict__ and dataset_info.sc_assess and os.path.exists(os.path.join(output_dir, "scaffolds.fasta")) and latest_found:
        latest_ctg = os.path.join(contig_dir, "latest_scaffolds.fasta")
        if not os.path.exists(latest_ctg):
            import glob
            prev_contigs = sorted(glob.glob(os.path.join(contig_dir, "*_scafs.fasta")))

            if len(prev_contigs) == 0:
                print("No previous scaffolds detected")
                latest_found = False
            else:
                os.chdir(contig_dir)
                os.symlink(os.path.basename(prev_contigs[-1]), "latest_scaffolds.fasta")
                os.chdir(working_dir)
        
        if latest_found:
            contigs = "scaffolds.fasta"
            quast_output_dir = os.path.join(output_dir, "QUAST_RESULTS_CMP_SC")
            qcode = run_quast(dataset_info, [os.path.join(output_dir, contigs), latest_ctg], quast_output_dir)

            #compare_misassemblies
            log.log("======= SCAFFOLD COMPARISON =======")
            if not cmp_misassemblies(quast_output_dir, "latest_scaffolds", os.path.splitext(contigs)[0]):
                rewrite_latest = False
                #exit_code = 14

    #writing log
    write_log(history_log, new_log, output_dir, dataset_info)


    #saving contigs
    if contig_dir != '':
        quast_contig_dir = os.path.join(contig_dir, "quast_results")
        if not os.path.exists(quast_contig_dir):
            os.makedirs(quast_contig_dir)

        name_prefix = datetime.datetime.now().strftime('%Y%m%d-%H%M')
#        if len(sys.argv) == 3:
#            name_prefix += "_" + sys.argv[2]
        print("Contigs have prefix " + name_prefix)

        if 'dipspades' in dataset_info.__dict__ and dataset_info.dipspades:
            shutil.copy(os.path.join(output_dir, "consensus_contigs.fasta"), os.path.join(contig_dir, name_prefix + "_ctg.fasta"))
        elif 'truseq' in dataset_info.__dict__ and dataset_info.truseq:
            shutil.copy(os.path.join(output_dir, "truseq_long_reads.fasta"), os.path.join(contig_dir, name_prefix + "_ctg.fasta"))
        else:
            shutil.copy(os.path.join(output_dir, "contigs.fasta"), os.path.join(contig_dir, name_prefix + "_ctg.fasta"))
        print("Contigs saved to " + os.path.join(contig_dir, name_prefix + "_ctg.fasta"))

        if rewrite_latest:
            print("Creating latest symlink")
            os.chdir(contig_dir)
            if os.path.islink("latest_contigs.fasta"):
                os.remove("latest_contigs.fasta")
            os.symlink(name_prefix + "_ctg.fasta", "latest_contigs.fasta")
            os.chdir(working_dir)

        scafs = os.path.join(output_dir, "scaffolds.fasta")
        if os.path.exists(scafs):
            shutil.copy(scafs, os.path.join(contig_dir, name_prefix + "_scafs.fasta"))
            print("Scaffolds saved to " + os.path.join(contig_dir, name_prefix + "_scafs.fasta"))

            if rewrite_latest:
                print("Creating latest symlink")
                os.chdir(contig_dir)
                if os.path.islink("latest_scaffolds.fasta"):
                    os.remove("latest_scaffolds.fasta")
                os.symlink(name_prefix + "_scafs.fasta", "latest_scaffolds.fasta")
                os.chdir(working_dir)

        before_rr = os.path.join(output_dir, "before_rr.fasta")
        if os.path.exists(before_rr):
            shutil.copy(before_rr, os.path.join(contig_dir, name_prefix + "_before_rr.fasta"))
            print("Contigs before resolve saved to " + os.path.join(contig_dir, name_prefix + "_before_rr.fasta"))

        preliminary_ctgs = os.path.join(output_dir, "first_pe_contigs.fasta")
        if os.path.exists(preliminary_ctgs):
            shutil.copy(before_rr, os.path.join(contig_dir, name_prefix + "_prelim.fasta"))
            print("Preliminary contigs after first RR saved to " + os.path.join(contig_dir, name_prefix + "_prelim.fasta"))

    #    before_corr = os.path.join(output_dir, "assembled_contigs.fasta")
    #    if os.path.exists(before_corr):
    #        shutil.copy(before_corr, os.path.join(contig_dir, name_prefix + "_before_corr.fasta"))

    #    scaff_before_corr = os.path.join(output_dir, "assembled_scaffolds.fasta")
    #    if os.path.exists(scaff_before_corr):
    #        shutil.copy(scaff_before_corr, os.path.join(contig_dir, name_prefix + "_scaff_before_corr.fasta"))

    #    boken_scaff = os.path.join(output_dir, "broken_scaffolds.fasta")
    #    if os.path.exists(boken_scaff):
    #        shutil.copy(boken_scaff, os.path.join(contig_dir, name_prefix + "_scaff_broken.fasta"))


    #    if quast_cmd != "":
    #        import glob
    #        sys.path.append('./src/tools/contig_analysis/')
    #        import compare_fasta
    #        print("Running quast for all saved contigs now...")
    #        os.system(quast_cmd + " -o " + quast_contig_dir + " " + os.path.join(contig_dir, "*.fasta") + " > " + os.path.join(contig_dir, "quast.log") + " 2> " + os.path.join(contig_dir, "quast.err"))
    #        print("Done")
    log.print_log()
    os.system("chmod -R 777 " + output_dir)
    sys.exit(exit_code)

except BaseException as e:
    print_exc()
    sys.exit(239)
