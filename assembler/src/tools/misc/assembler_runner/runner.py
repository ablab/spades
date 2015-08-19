#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import shutil
import subprocess
import sys
import getopt
import logging
from datetime import datetime

# binaries
runner_dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
quast_bin = 'quast' # '/home/alex/biolab/master_rep/quast/quast.py'
abyss_base = '/Molly/oldstorage/acestorage/software/abyss_1.3.6/abyss_build/bin/'
msrca_base = '/Molly/oldstorage/acestorage/software/MaSuRCA-2.0.3.1/bin/'
msrca_config = os.path.join(runner_dirpath, 'MSRCA_config.txt')
msrca_mp_config = os.path.join(runner_dirpath, 'MSRCA_mp_config.txt')
soap_base = '/Molly/oldstorage/acestorage/software/SOAPdenovo_2.04/'
soap_run_script = os.path.join(runner_dirpath, 'SOAP_assemble.sh')
soap_config = os.path.join(runner_dirpath, 'SOAP_config.txt')
soap_mp_config = os.path.join(runner_dirpath, 'SOAP_mp_config.txt')
idba_base = '/Molly/oldstorage/acestorage/software/idba-1.1.1/bin/'
spades_base = '/Molly/oldstorage/acestorage/software/SPAdes-3.5.0-Linux/bin/'
ray_base = '/Molly/oldstorage/acestorage/software/Ray-v2.2.0/Ray-Large-k-mers/'
cabog_base = '/Molly/oldstorage/acestorage/software/CABOG_wgs_7.0/Linux-amd64/bin'
cabog_run_script = os.path.join(runner_dirpath, 'CABOG_assemble.sh')
cabog_mp_run_script = os.path.join(runner_dirpath, 'CABOG_mp_assemble.sh')
cabog_config = os.path.join(runner_dirpath, 'CABOG_config.txt')
sga_base = '/Molly/oldstorage/acestorage/software/SGA-0.10.10/build/bin/'
sga_run_script = os.path.join(runner_dirpath, 'SGA_assemble.sh')
sga_mp_run_script = os.path.join(runner_dirpath, 'SGA_mp_assemble.sh')
sga_cor_run_script = os.path.join(runner_dirpath, 'SGA_assemble_cor.sh')
sga_mp_cor_run_script = os.path.join(runner_dirpath, 'SGA_mp_assemble_cor.sh')
velvet_base = '/Molly/oldstorage/acestorage/software/velvet_1.2.10/'
velvet_run_script = os.path.join(runner_dirpath, 'VELVET_assemble.sh')
velvet_mp_run_script = os.path.join(runner_dirpath, 'VELVET_mp_assemble.sh')
velvet_sc_base = '/Molly/oldstorage/acestorage/software/velvet-sc/'
velvet_sc_run_script = os.path.join(runner_dirpath, 'VELVETSC_assemble.sh')
velvet_sc_mp_run_script = os.path.join(runner_dirpath, 'VELVETSC_mp_assemble.sh')
mira_base = '/Molly/oldstorage/acestorage/software/mira-4.0rc1/build/bin/'
mira_config = os.path.join(runner_dirpath, 'MIRA_config.txt')
mira_mp_config = os.path.join(runner_dirpath, 'MIRA_mp_config.txt')
perga_base = '/Molly/oldstorage/acestorage/software/perga-0.5.03.01/PERGA/bin/'

# IMPORTANT CONSTANTS
long_options = "12= mp12= min-k= max-k= step-k= threads= is= dev= assemblers= mira-tmp-dir= corrected sc mp1= mp2= mp-orient= mp-is= mp-dev= phred64".split()
short_options = "o:1:2:r:t:"

# options
assemblers_to_run = ['abyss', 'msrca', 'soap', 'idba', 'spades', 'ray', 'cabog', 'sga', 'velvet', 'velvet-sc', 'mira', 'perga']
output_dir = None
mira_tmp_dir = os.path.abspath(os.path.expanduser('~/mira_tmp'))
left = None
right = None
interlaced = None
mate_left = None
mate_right = None
mate_interlaced = None
mate_orient = 'rf'
mate_is = None
mate_is_dev = None
insert_size = None
is_dev = None
max_threads = 8
reference = None
corrected = False
single_cell = False
phred64 = False

k_2level_selection = True
min_k_rough = 35
max_k_rough = 95
step_k_rough = 20
delta_k_accurate = 10
step_k_accurate = 4

# aux directories prefixes
BEST_DIR = "best"
ARCHIVE_DIR = "archive"
# aux suffixes
CONTIGS_SUFF = ".ctg.fasta"
SCAFFOLDS_SUFF = ".scf.fasta"
LOG_SUFF = ".log.txt"
QUAST_REPORT_SUFF = ".report.txt"

# for assemblers with configs
PARAM_PREFIX = "RUNNER_PARAM_"
common_params_subst_dict = dict()

# for MIRA (should crop long read names)
MAX_READ_NAME_LENGTH = 30
ALLOWED_READ_NAME_PREFIX = 8


def error(err_str, prefix="ERROR: ", stdout=False):
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    prefix = current_time + " " + prefix
    if stdout:
        sys.stderr.write("\n" + prefix + " " + err_str + "\n\n")
        sys.stderr.flush()
        sys.exit(1)
    else:
        log = logging.getLogger('runner')
        log.info("\n" + prefix + " " + err_str + "\n")


def warning(warn_str, prefix="WARNING: ", stdout=False):
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    prefix = current_time + " " + prefix
    if stdout:
        sys.stdout.write("\n" + prefix + " " + warn_str + "\n\n")
        sys.stdout.flush()
    else:
        log = logging.getLogger('runner')
        log.info("\n" + prefix + " " + warn_str + "\n")


def info(info_str, prefix=""):
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    prefix = current_time + " " + prefix
    log = logging.getLogger('runner')
    log.info(prefix + " " + info_str)


def get_max_id(values):
    max_val = max(values)
    for id, v in enumerate(values):
        if v == max_val:
            return id
    return None


def get_int(value):
    try:
        return int(value)
    except ValueError:
        return -1


def check_options():
    global output_dir
    global left
    global right
    global interlaced
    global reference

    if not output_dir:
        error("output_dir is not set!", stdout=True)
        sys.exit(1)
    if (not left or not right) and not (interlaced and assemblers_to_run[0] == assemblers_to_run[-1] == 'velvet-sc'):
        error("you should set both left and right reads!", stdout=True)
        sys.exit(1)
    if not insert_size or not is_dev:
        error("you should set both insert size mean (--is) and insert size deviation (--dev)!", stdout=True)
        sys.exit(1)
    if interlaced and assemblers_to_run[0] == assemblers_to_run[-1] == 'velvet-sc':
        if not os.path.isfile(interlaced):
            error("file with interlaced reads doesn't exist! " + interlaced, stdout=True)
            sys.exit(1)
    else:
        if not os.path.isfile(left):
            error("file with left reads doesn't exist! " + left, stdout=True)
            sys.exit(1)
        if not os.path.isfile(right):
            error("file with right reads doesn't exist! " + right, stdout=True)
            sys.exit(1)
    if not reference or not os.path.isfile(reference):
        warning("reference is NOT SET, Quast will detect best assembly by N50!", stdout=True)
    else:
        reference = os.path.abspath(reference)
    output_dir = os.path.abspath(output_dir)
    if left and right:
        left = os.path.abspath(left)
        right = os.path.abspath(right)
    else:
        left = ''
        right = ''
    if interlaced:
        interlaced = os.path.abspath(interlaced)


def fill_in_common_params_subst_dict():
    global common_params_subst_dict
    common_params_subst_dict['INSERT_SIZE'] = str(insert_size)
    common_params_subst_dict['DEVIATION'] = str(is_dev)
    common_params_subst_dict['LEFT'] = left
    common_params_subst_dict['RIGHT'] = right
    common_params_subst_dict['THREADS'] = str(max_threads)


def update_runner_params(filename, params_subst_dict):
    old_file = open(filename, 'r')
    old_lines = old_file.readlines()
    old_file.close()
    new_file = open(filename, 'w')
    for line in old_lines:
        for k, v in params_subst_dict.items():
            to_check = PARAM_PREFIX + k.upper()
            if to_check in line:
                line = line.replace(to_check, v)
        new_file.write(line)
    new_file.close()


# run assemblers functions
def run_K_dependent_assembler(name, K, run_tool_function, err_number):
    info("\tstart " + name + " with K=" + str(K))
    out_dir = os.path.join(output_dir, name, "K" + str(K))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    archive_dir = os.path.join(output_dir, ARCHIVE_DIR)
    cwd = os.getcwd()
    os.chdir(out_dir)
    log_filename = os.path.join(out_dir, name + LOG_SUFF)
    return_code, contigs_path, scaffolds_path = run_tool_function(K, log_filename)
    if return_code:
        warning(name + " with K=" + str(K) + " returned non-zero code!")

    result = None
    for (type, src, suff) in [("contigs", contigs_path, CONTIGS_SUFF), ("scaffolds", scaffolds_path, SCAFFOLDS_SUFF)]:
        if src and os.path.isfile(src):
            result = os.path.join(archive_dir, name + ".K" + str(K) + suff)
            shutil.copy(src, result)
        else:
            if src:
                warning(type + " not found! " + os.path.abspath(src))
            else:
                warning(type + " not found! (assembler returned None)")
    shutil.copy(log_filename, os.path.join(archive_dir, name + ".K" + str(K) + LOG_SUFF))
    if result or (err_number > 0):  #TODISCUSS: maybe not 0 here
        shutil.rmtree(out_dir)
    os.chdir(cwd)
    info("\tfinish " + name + " with K=" + str(K))
    return result


# run tools for assemblers depended on K
def run_abyss(K, log_filename):
    log_file = open(log_filename, 'a')
    cmd_line = os.path.join(abyss_base, "abyss-pe") + " k=%d n=5 s=100 name=asm lib='reads' reads='%s %s'" % \
        (K, left, right) + " -j " + str(max_threads)
    if mate_left and mate_right:
        cmd_line += " mp_lib='mp_reads' mp_reads='%s %s' " % (mate_left, mate_right)
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.flush()
    #cmd_line = '/home/alex/biolab/TEST/RUNNER/data/abyss/abyss_sim.sh %s %d %s' % ('abyss', K, os.path.join(output_dir, 'abyss', "K" + str(K)))
    #cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    return_code = os.system(cmd_line)
    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "asm-contigs.fa"
    scaffolds_path = "asm-scaffolds.fa"
    return return_code, contigs_path, scaffolds_path


def run_msrca(K, log_filename):
    log_file = open(log_filename, 'a')
    cur_msrca_config = os.path.join(os.getcwd(), 'config.txt')
    #shutil.copy(msrca_config, cur_msrca_config)
    params_subst_dict = dict(common_params_subst_dict)
    params_subst_dict['K'] = str(K)
    if mate_left and mate_right and mate_is and mate_is_dev:
        params_subst_dict['MP_LEFT'] = mate_left
        params_subst_dict['MP_RIGHT'] = mate_right
        params_subst_dict['MP_INSERT_SIZE'] = str(mate_is)
        if mate_orient == 'rf':
            params_subst_dict['MP_DEVIATION'] = str(mate_is_dev)
        else:
            params_subst_dict['MP_DEVIATION'] = str(-mate_is_dev)
        shutil.copy(msrca_mp_config, cur_msrca_config)
    else:
        shutil.copy(msrca_config, cur_msrca_config)
    update_runner_params(cur_msrca_config, params_subst_dict)
    cmd_line = os.path.join(msrca_base, "runSRCA.pl") + " " + cur_msrca_config
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started step 1: " + cmd_line + "\n\n")
    log_file.write("Config content:\n\n")
    for line in open(cur_msrca_config, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    if return_code:
        warning("MSRCA generating step returned non-zero code!")
    assemble_script = os.path.abspath(os.path.join(os.getcwd(), 'assemble.sh'))
    if not os.path.isfile(assemble_script):
        warning("MSRCA generating step didn't create assemble.sh!")
        return 0, None, None

    cmd_line = assemble_script
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started step 2: " + cmd_line + "\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    ###
    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "CA/10-gapclose/genome.ctg.fasta"
    scaffolds_path = "CA/10-gapclose/genome.scf.fasta"
    return return_code, contigs_path, scaffolds_path


def run_soap(K, log_filename):
    log_file = open(log_filename, 'a')
    cur_soap_config = os.path.abspath(os.path.join(os.getcwd(), 'config.txt'))
    #shutil.copy(soap_config, cur_soap_config)
    cur_soap_run_script = os.path.abspath(os.path.join(os.getcwd(), 'run.sh'))
    shutil.copy(soap_run_script, cur_soap_run_script)
    params_subst_dict = dict(common_params_subst_dict)
    params_subst_dict['K'] = str(K)
    params_subst_dict['SOAP_BASE'] = soap_base
    params_subst_dict['SOAP_CONFIG'] = cur_soap_config
    if mate_left and mate_right and mate_is:
        params_subst_dict['MP_LEFT'] = mate_left
        params_subst_dict['MP_RIGHT'] = mate_right
        params_subst_dict['MP_INSERT_SIZE'] = str(mate_is)
        shutil.copy(soap_mp_config, cur_soap_config)
    else:
        shutil.copy(soap_config, cur_soap_config)
    update_runner_params(cur_soap_config, params_subst_dict)
    update_runner_params(cur_soap_run_script, params_subst_dict)
    cmd_line = cur_soap_run_script
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.write("Run script content:\n")
    for line in open(cur_soap_run_script, 'r'):
        log_file.write(line)
    log_file.write("\n\nConfig content:\n")
    for line in open(cur_soap_config, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()

    return_code = os.system(cmd_line)

    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "asm.contig"
    scaffolds_path = "asm.new.scafSeq"
    return return_code, contigs_path, scaffolds_path


def run_velvet(K, log_filename):
    log_file = open(log_filename, 'a')
    cur_velvet_run_script = os.path.abspath(os.path.join(os.getcwd(), 'run.sh'))
    #shutil.copy(velvet_run_script, cur_velvet_run_script)
    params_subst_dict = dict(common_params_subst_dict)
    params_subst_dict['K'] = str(K)
    params_subst_dict['VELVET_BASE'] = velvet_base
    params_subst_dict['OUTPUT_DIR'] = os.path.dirname(log_filename)
    if mate_left and mate_right and mate_is and mate_is_dev:
        params_subst_dict['MP_LEFT'] = mate_left
        params_subst_dict['MP_RIGHT'] = mate_right
        params_subst_dict['MP_INSERT_SIZE'] = str(mate_is)
        params_subst_dict['MP_DEVIATION'] = str(mate_is_dev)
        shutil.copy(velvet_mp_run_script, cur_velvet_run_script)
    else:
        shutil.copy(velvet_run_script, cur_velvet_run_script)
    update_runner_params(cur_velvet_run_script, params_subst_dict)
    cmd_line = cur_velvet_run_script
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.write("Run script content:\n")
    for line in open(cur_velvet_run_script, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()

    return_code = os.system(cmd_line)

    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "contigs.fa"
    scaffolds_path = "scaffolds.fa"
    return return_code, contigs_path, scaffolds_path


def run_velvet_sc(K, log_filename):
    log_file = open(log_filename, 'a')
    cur_velvet_run_script = os.path.abspath(os.path.join(os.getcwd(), 'run.sh'))
    #shutil.copy(velvet_run_script, cur_velvet_run_script)
    params_subst_dict = dict(common_params_subst_dict)
    params_subst_dict['K'] = str(K)
    params_subst_dict['VELVET_BASE'] = velvet_sc_base
    params_subst_dict['OUTPUT_DIR'] = os.path.dirname(log_filename)
    params_subst_dict['INTERLACED'] = interlaced if interlaced else ''
    if mate_left and mate_right and mate_is and mate_is_dev:
        params_subst_dict['MP_LEFT'] = mate_left
        params_subst_dict['MP_RIGHT'] = mate_right
        params_subst_dict['MP_INSERT_SIZE'] = str(mate_is)
        params_subst_dict['MP_DEVIATION'] = str(mate_is_dev)
        params_subst_dict['MP_INTERLACED'] = mate_interlaced if mate_interlaced else ''
        shutil.copy(velvet_sc_mp_run_script, cur_velvet_run_script)
    else:
        shutil.copy(velvet_sc_run_script, cur_velvet_run_script)
    update_runner_params(cur_velvet_run_script, params_subst_dict)
    cmd_line = cur_velvet_run_script
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.write("Run script content:\n")
    for line in open(cur_velvet_run_script, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()

    return_code = os.system(cmd_line)

    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "contigs.fa"
    scaffolds_path = "scaffolds.fa"
    return return_code, contigs_path, scaffolds_path


def run_sga(K, log_filename):
    log_file = open(log_filename, 'a')
    cur_sga_run_script = os.path.abspath(os.path.join(os.getcwd(), 'run.sh'))
    if corrected:
        if mate_left and mate_right:
            shutil.copy(sga_mp_cor_run_script, cur_sga_run_script)
        else:
            shutil.copy(sga_cor_run_script, cur_sga_run_script)
    else:
        if mate_left and mate_right:
            shutil.copy(sga_mp_run_script, cur_sga_run_script)
        else:
            shutil.copy(sga_run_script, cur_sga_run_script)
    params_subst_dict = dict(common_params_subst_dict)
    params_subst_dict['K'] = str(K)
    params_subst_dict['SGA_BASE'] = sga_base
    update_runner_params(cur_sga_run_script, params_subst_dict)
    cmd_line = cur_sga_run_script
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    target_link = 'frag1'
    log_file.write("Started preparation 1: creating symlink from %s to %s\n\n" % (left, target_link))
    os.symlink(left, target_link)
    target_link = 'frag2'
    log_file.write("Started preparation 2: creating symlink from %s to %s\n\n" % (right, target_link))
    os.symlink(right, target_link)
    if mate_left and mate_right:
        target_link = 'jump1'
        log_file.write("Started preparation 3: creating symlink from %s to %s\n\n" % (mate_left, target_link))
        os.symlink(mate_left, target_link)
        target_link = 'jump2'
        log_file.write("Started preparation 4: creating symlink from %s to %s\n\n" % (mate_right, target_link))
        os.symlink(mate_right, target_link)
    log_file.write("Started assembling: " + cmd_line + "\n\n")
    log_file.write("Run script content:\n")
    for line in open(cur_sga_run_script, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()

    return_code = os.system(cmd_line)

    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "ctg.fasta"
    scaffolds_path = "scf.fasta"
    return return_code, contigs_path, scaffolds_path


def run_ray(K, log_filename):
    log_file = open(log_filename, 'a')
    out_subdir = "out"
    if os.path.isdir(out_subdir):
        shutil.rmtree(out_subdir)
    cmd_line = "mpirun -np %d" % (max_threads / 2) + " " + os.path.join(ray_base, "Ray") + " -o out -k %d -p %s %s" % \
                                                      (K, left, right)  ## TODO max_threads / 2 or not.
    if mate_left and mate_right:
        if mate_orient == 'fr':  # TODO: else
            cmd_line += " -p %s %s " % (mate_left, mate_right)
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = os.path.join(out_subdir, "Contigs.fasta")
    scaffolds_path = os.path.join(out_subdir, "Scaffolds.fasta")
    return return_code, contigs_path, scaffolds_path


def run_perga(K, log_filename):
    log_file = open(log_filename, 'a')
    cmd_line = os.path.join(perga_base, "perga") + " all -k %d -p 1 -f %s %s -d . -ins_len %d -ins_sdev %d " % \
        (K, left, right, insert_size, is_dev)
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "contigs.fa"
    scaffolds_path = "scaffolds.fa"
    return return_code, contigs_path, scaffolds_path


def iterate_K(name, run_tool_function):
    info("START " + name)
    out_dir = os.path.join(output_dir, name)
    archive_dir = os.path.join(output_dir, ARCHIVE_DIR)
    best_dir = os.path.join(output_dir, BEST_DIR)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    global_err_number = 0
    best_K = None
    best_assembly = None
    for K_level_selection in ([1, 2] if k_2level_selection else [1]):
        # START OF TOOL with different K
        if K_level_selection == 1:
            k_list = range(min_k_rough, max_k_rough + 1, step_k_rough)
        else:
            k_2level_seed = best_K # best K from first level
            k_list = sorted([k_2level_seed - delta_k_accurate, k_2level_seed - delta_k_accurate + step_k_accurate,
                             k_2level_seed + delta_k_accurate - step_k_accurate, k_2level_seed + delta_k_accurate])
            while best_K in k_list:
                k_list.remove(best_K)
        err_number = 0
        results_dict = dict()
        if K_level_selection == 2:
            results_dict[best_K] = best_assembly

        for K in k_list:
            result = run_K_dependent_assembler(name, K, run_tool_function, err_number)
            if not result:
                err_number += 1
            else:
                results_dict[K] = result

        # START OF QUAST evaluating
        quast_out_dir = os.path.join(out_dir, "quast_results")
        cmd_line = quast_bin + " " + " ".join(results_dict.values()) + " -o " + quast_out_dir
        key_metric = "N50"
        aux_metric = "N50"
        if reference:
            cmd_line += " -R " + reference
            key_metric = "NGA50"
        #info("running QUAST: " + cmd_line)
        info("running QUAST (%d assemblies)" % len(results_dict.values()))
        subprocess.call(cmd_line.split(), stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
        quast_report = os.path.join(quast_out_dir, "report.txt")
        if os.path.isfile(quast_report):
            shutil.copy(quast_report, os.path.join(archive_dir, name + ".K_lev_" + str(K_level_selection) + QUAST_REPORT_SUFF))
            assemblies_ids = []
            assemblies_results = []
            assemblies_aux_results = []
            for line in open(quast_report, 'r'):
                if line.startswith("Assembly"):
                    assemblies_ids = line.split()[1:]
                if line.startswith(key_metric):
                    assemblies_results = map(get_int, line.split()[1:])
                if line.startswith(aux_metric):
                    assemblies_aux_results = map(get_int, line.split()[1:])
            if len(assemblies_ids) != len(results_dict.keys()):
                warning("QUAST did NOT report results for ALL assemblies!")
            if not assemblies_results:
                if assemblies_aux_results:
                    warning("QUAST did NOT report results for metric " + key_metric + ", using " + aux_metric)
                    assemblies_results = assemblies_aux_results
                    key_metric = aux_metric
                else:
                    warning("QUAST failed to evaluates assemblies, ABORTING with " + name + ", see details in: " + quast_out_dir)
                    return
            best_assembly = assemblies_ids[get_max_id(assemblies_results)]
            best_K = int(best_assembly.split('.')[-2][1:]) # [-3] is "Kxx" fragment
            info("best K is " + str(best_K) + " (best assembly is " + best_assembly +
                 ", its result (" + key_metric + ") is " + str(max(assemblies_results)) + ")")
        else:
            warning("QUAST failed to generate report, ABORTING with " + name + ", see details in: " + quast_out_dir)
            return
        # END OF QUAST

        if best_K in results_dict.keys():
            best_assembly = results_dict[best_K]
            if K_level_selection == 2 or not k_2level_selection:
                shutil.copy(best_assembly, best_dir)
                info("best assembly copied to " + best_dir)
                if best_assembly.endswith(SCAFFOLDS_SUFF):
                    best_contigs = best_assembly[:len(best_assembly) - len(SCAFFOLDS_SUFF)] + CONTIGS_SUFF
                    if os.path.isfile(best_contigs):
                        shutil.copy(best_contigs, best_dir)
            else:
                info("best K selected for the 2nd level of K iterating: " + str(best_K))
        else:
            warning("something went wrong with detecting best assembly, best K is " + str(best_K) +
                    ", but results_dict contains results only for " + str(results_dict.keys()) + "! ABORTING with " + name)
            return
        global_err_number += err_number

    if not global_err_number:
        shutil.rmtree(out_dir)
    info("FINISH " + name + "\n")


def run_K_independent_assembler(name, run_tool_function):
    info("START " + name)
    out_dir = os.path.join(output_dir, name)
    archive_dir = os.path.join(output_dir, ARCHIVE_DIR)
    best_dir = os.path.join(output_dir, BEST_DIR)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    cwd = os.getcwd()
    os.chdir(out_dir)
    log_filename = os.path.join(out_dir, name + LOG_SUFF)
    if name == 'sga': # in corrected mode, SGA is K independent de facto
        return_code, contigs_path, scaffolds_path = run_tool_function(0, log_filename)
    else:
        return_code, contigs_path, scaffolds_path = run_tool_function(log_filename)
    if return_code:
        warning(name + " returned non-zero code!")
    result = None
    for (type, src, suff) in [("contigs", contigs_path, CONTIGS_SUFF), ("scaffolds", scaffolds_path, SCAFFOLDS_SUFF)]:
        if src and os.path.isfile(src):
            result = os.path.join(archive_dir, name + suff)
            shutil.copy(src, result)
        else:
            if src:
                warning(type + " not found! " + os.path.abspath(src))
            else:
                warning(type + " not found! (assembler returned None)")
    shutil.copy(log_filename, os.path.join(archive_dir, name + LOG_SUFF))
    os.chdir(cwd)

    if result and os.path.isfile(result):
        shutil.rmtree(out_dir)
        best_assembly = result
        shutil.copy(best_assembly, best_dir)
        info("best assembly copied to " + best_dir)
        if best_assembly.endswith(SCAFFOLDS_SUFF):
            best_contigs = best_assembly[:len(best_assembly) - len(SCAFFOLDS_SUFF)] + CONTIGS_SUFF
            if os.path.isfile(best_contigs):
                shutil.copy(best_contigs, best_dir)
    else:
        warning("something went wrong and " + name + " didn't generate neither contigs nor scaffolds! ABORTING with " + name)
        return
    info("FINISH " + name + "\n")


# run tools for assemblers independed on K
def run_spades(log_filename):
    log_file = open(log_filename, 'a')
    cmd_line = os.path.join(spades_base, "spades.py") + " --careful -o . -1 %s -2 %s -t %d" % (left, right, max_threads)
    if corrected:
        cmd_line += " --only-assembler"
    if single_cell:
        cmd_line += " --sc"
    if mate_left and mate_right:
        cmd_line += " --mp1-1 %s --mp1-2 %s --mp1-%s " % (mate_left, mate_right, mate_orient)
    #cmd_line = os.path.join(spades_base, "spades.py") + " --sc -o . -1 %s -2 %s -t %d" % (left, right, max_threads)
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "contigs.fasta"
    scaffolds_path = "scaffolds.fasta"
    return return_code, contigs_path, scaffolds_path


def run_idba(log_filename):
    log_file = open(log_filename, 'a')
    shuffled_reads = os.path.abspath('shuffled.fasta')
    cmd_line = os.path.join(idba_base, "fq2fa") + " --merge %s %s %s" % (left, right, shuffled_reads)
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started preparing step: " + cmd_line + "\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    if return_code:
        warning("IDBA preparing step (shuffling two fastQ into one fastA) returned non-zero code!")
    if not os.path.isfile(shuffled_reads):
        warning("IDBA preparing step didn't create shuffled reads!")
        return 0, None, None

    shuffled_mp_reads = None
    if mate_left and mate_right: 
        if mate_orient == 'fr':  # TODO: else
            shuffled_mp_reads = os.path.abspath('shuffled_mp.fasta')
            cmd_line = os.path.join(idba_base, "fq2fa") + " --merge %s %s %s" % (mate_left, mate_right, shuffled_mp_reads)
            cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
            log_file.write("Started preparing step-2 (mate-pairs): " + cmd_line + "\n\n")
            log_file.flush()
            return_code = os.system(cmd_line)
            if return_code:
                warning("IDBA preparing step-2 (shuffling two fastQ into one fastA) returned non-zero code!")
            if not os.path.isfile(shuffled_mp_reads):
                warning("IDBA preparing step-2 didn't create shuffled reads!")
                return 0, None, None

    cmd_line = os.path.join(idba_base, "idba_ud") + " -o out -r %s --num_threads %d" % (shuffled_reads, max_threads)
    if shuffled_mp_reads:
        cmd_line += " --read_level_2 %s " % (shuffled_mp_reads)
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)
    log_file.write("Started assembling step: " + cmd_line + "\n\n")
    log_file.flush()
    return_code = os.system(cmd_line)
    ###
    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "out/contig.fa"
    scaffolds_path = "out/scaffold.fa"
    if shuffled_mp_reads:
        scaffolds_path = "out/scaffold-level-2.fa"
    return return_code, contigs_path, scaffolds_path


def run_cabog(log_filename):
    log_file = open(log_filename, 'a')
    cur_cabog_config = os.path.abspath(os.path.join(os.getcwd(), 'config.txt'))
    shutil.copy(cabog_config, cur_cabog_config)
    cur_cabog_run_script = os.path.abspath(os.path.join(os.getcwd(), 'run.sh'))
    #shutil.copy(cabog_run_script, cur_cabog_run_script)
    params_subst_dict = dict(common_params_subst_dict)
    params_subst_dict['CABOG_BASE'] = cabog_base
    params_subst_dict['CABOG_CONFIG'] = cur_cabog_config
    params_subst_dict['STDOUT_LOG'] = log_filename
    params_subst_dict['STDERR_LOG'] = log_filename
    if phred64:
        params_subst_dict['PHRED_TYPE'] = 'illumina'
    else:
        params_subst_dict['PHRED_TYPE'] = 'sanger'
    if mate_left and mate_right and mate_is and mate_is_dev:
        log_file.write("Starting with mate-pairs\n")
        params_subst_dict['MP_LEFT'] = mate_left
        params_subst_dict['MP_RIGHT'] = mate_right
        params_subst_dict['MP_INSERT_SIZE'] = str(mate_is)
        params_subst_dict['MP_DEVIATION'] = str(mate_is_dev)
        if mate_orient == 'fr':
            params_subst_dict['MP_ORIENT'] = '-innie'
        else:
            params_subst_dict['MP_ORIENT'] = '-outtie'
        shutil.copy(cabog_mp_run_script, cur_cabog_run_script)
    else:
        log_file.write("Starting only with paired-end\n")
        shutil.copy(cabog_run_script, cur_cabog_run_script)
    update_runner_params(cur_cabog_run_script, params_subst_dict)
    cmd_line = cur_cabog_run_script
    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.write("Run script content:\n")
    for line in open(cur_cabog_run_script, 'r'):
        log_file.write(line)
    log_file.write("\n\nConfig content:\n")
    for line in open(cur_cabog_config, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()

    return_code = os.system(cmd_line)

    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    contigs_path = "9-terminator/asm.ctg.fasta"
    scaffolds_path = "9-terminator/asm.scf.fasta"
    return return_code, contigs_path, scaffolds_path


def run_mira(log_filename):
    def crop_read_names_if_needed(left_filename, right_filename, dir_name_to_save):
        crop_is_needed = False
        for id, line in enumerate(open(left_filename, 'r')):  # names in right_filename should be the same!
            if (id % 4) == 0:
                read_name = line.split()[0]
                if len(read_name) > MAX_READ_NAME_LENGTH or not read_name.endswith('/1'):
                    crop_is_needed = True
                    break
        if not crop_is_needed:
            return left_filename, right_filename
        else:
            new_left_filename = os.path.join(dir_name_to_save, os.path.basename(left_filename))
            new_right_filename = os.path.join(dir_name_to_save, os.path.basename(right_filename))
            for (new_fname, old_fname, suffix) in [(new_left_filename, left_filename, '/1'),
                (new_right_filename, right_filename, '/2')]:
                new_file = open(new_fname, 'w')
                for id, line in enumerate(open(old_fname, 'r')):
                    if (id % 4) == 0:
                        new_read_name = line.split()[0]
                        if not new_read_name.endswith(suffix):
                            new_read_name += suffix
                        if len(new_read_name) > MAX_READ_NAME_LENGTH:
                            new_read_name = new_read_name[:1 + ALLOWED_READ_NAME_PREFIX] + \
                                new_read_name[len(new_read_name) - (MAX_READ_NAME_LENGTH - ALLOWED_READ_NAME_PREFIX):]
                        new_read_name += " " + " ".join(line.split()[1:])
                        new_file.write(new_read_name + '\n')
                    else:
                        new_file.write(line)
                new_file.close()
            return new_left_filename, new_right_filename

    log_file = open(log_filename, 'a')
    cur_mira_config = os.path.abspath(os.path.join(os.getcwd(), 'config.txt'))
    params_subst_dict = dict(common_params_subst_dict)
    
    #shutil.copy(mira_config, cur_mira_config)
    cropped_left, cropped_right = crop_read_names_if_needed(left, right, os.path.dirname(log_filename))
    if mate_left and mate_right and mate_is and mate_is_dev:
        log_file.write("Starting with mate-pairs\n")
        cropped_mate_left, cropped_mate_right = crop_read_names_if_needed(mate_left, mate_right, os.path.dirname(log_filename))
        shutil.copy(mira_mp_config, cur_mira_config)
        params_subst_dict['MIRA_MP_MIN'] = str(max(0, mate_is - mate_is_dev)) # TODO: check this
        params_subst_dict['MIRA_MP_MAX'] = str(min(mate_is + mate_is_dev, 2 * mate_is)) # TODO: check this
        params_subst_dict['MATE_LEFT'] = cropped_mate_left
        params_subst_dict['MATE_RIGHT'] = cropped_mate_right
        params_subst_dict['MP_ORIENT'] = mate_orient.upper()
    else:
        log_file.write("Starting only with paired-ends\n")
        shutil.copy(mira_config, cur_mira_config)
    
    global mira_tmp_dir
    original_mira_tmp_dir = mira_tmp_dir
    i = 1
    while os.path.isdir(mira_tmp_dir):
        mira_tmp_dir = original_mira_tmp_dir + "_" + str(i)
        i += 1
    os.makedirs(mira_tmp_dir)
    params_subst_dict['MIRA_TMP_DIR'] = mira_tmp_dir
    params_subst_dict['MIRA_MIN'] = str(max(0, insert_size - 9 * is_dev / 2)) # TODO: check this
    params_subst_dict['MIRA_MAX'] = str(min(insert_size + 9 * is_dev / 2, 2 * insert_size)) # TODO: check this
    params_subst_dict['LEFT'] = cropped_left
    params_subst_dict['RIGHT'] = cropped_right
    update_runner_params(cur_mira_config, params_subst_dict)
    cmd_line = os.path.join(mira_base, 'mira') + " " + cur_mira_config
    cmd_line += " >>%s 2>>%s" % (log_filename, log_filename)

    log_file.write("Started: " + cmd_line + "\n\n")
    log_file.write("\n\nConfig content:\n")
    for line in open(cur_mira_config, 'r'):
        log_file.write(line)
    log_file.write("\n\n")
    log_file.flush()

    return_code = os.system(cmd_line)

    log_file.write("Finished, log is %s\n" % log_filename)
    log_file.close()
    shutil.rmtree(mira_tmp_dir)
    contigs_path = "MyFirstAssembly_assembly/MyFirstAssembly_d_results/MyFirstAssembly_out.unpadded.fasta"
    scaffolds_path = None # mira doesn't have a scaffolder
    return return_code, contigs_path, scaffolds_path


def main():
    global output_dir
    global left
    global right
    global interlaced
    global reference
    global insert_size
    global is_dev
    global max_threads
    global min_k_rough
    global max_k_rough
    global step_k_rough
    global k_2level_selection
    global assemblers_to_run
    global mira_tmp_dir
    global corrected
    global single_cell
    global mate_left
    global mate_right
    global mate_interlaced
    global mate_orient
    global mate_is
    global mate_is_dev
    global phred64

    if len(sys.argv) == 1:
        print("Assemblers runner v.1.1")
        print("Usage: " + str(sys.argv[0]) + " [options]")
        print("")
        print("Required options:")
        print("-o\t<output_dir>\tbase output dir")
        print("-1\t<file>\t\tleft reads (in FASTQ)")
        print("-2\t<file>\t\tright reads (in FASTQ)")
        print("--is\t<int>\t\tinsert size mean")
        print("--dev\t<int>\t\tinsert size standard deviation")
        print("Mate-pairs options:")
        print("--mp-orient\t<rf,fr>\tmate-pairs orientation (RF -- default or FR)")
        print("--mp1\t<file>\t\tleft mate-pairs reads (in FASTQ)")
        print("--mp2\t<file>\t\tright mate-pairs reads (in FASTQ)")
        print("--mp-is\t<int>\t\tmate-pairs insert size mean")
        print("--mp-dev\t<int>\t\tmate-pairs insert size standard deviation")
        print("\nRecommended options:")
        print("-r\t\t<file>\treference (in FASTA)")
        print("--threads\t<int>\tmax threads number (default is %d)" % max_threads)
        print("\nAssembler specific options:")
        print("--sc\t\t\t\tdataset is single-cell [for SPAdes]")
        print("--12\t<file>\t\tinterlaced reads (in FASTQ) [for Velvet-SC]")
        print("--mp12\t<file>\t\tinterlaced mate-pairs reads (in FASTQ) [for Velvet-SC]")
        print("--mira-tmp-dir\t<dir>\ttemporary dir [for MIRA assembler] (default is %s)" % mira_tmp_dir)
        print("--corrected\t\t\treads are corrected (disable error correction [for SGA, SPAdes])")
        print("\nAdvanced options:")
        print("--min-k\t\t<int>\tmin K (default is auto from %s)" % str(min_k_rough - delta_k_accurate))
        print("--max-k\t\t<int>\tmax K (default is auto up to %s)" % str(max_k_rough + delta_k_accurate))
        print("--step-k\t<int>\tstep for K changing (default is auto, first %s and than %s)" % (str(step_k_rough), str(step_k_accurate)))
        print("--phred64\t\tdataset quality encoding is Phred+64 (default is Phred+33)")
        print("--assemblers\t<str,str>\tlist of assemblers to run (default is %s)" % assemblers_to_run)

        sys.exit(0)

    try:
        options, not_options = getopt.gnu_getopt(sys.argv, short_options, long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        sys.stderr.flush()
        sys.exit(1)

    for opt, arg in options:
        if opt == '-o':
            output_dir = arg
        elif opt == '-1':
            left = arg
        elif opt == '-2':
            right = arg
        elif opt == '--12':
            interlaced = arg
        elif opt == '-r':
            reference = arg
        elif opt == '--is':
            insert_size = int(arg)
        elif opt == '--dev':
            is_dev = int(arg)
        elif opt == '--mp1':
            mate_left = arg
        elif opt == '--mp2':
            mate_right = arg
        elif opt == '--mp12':
            mate_interlaced = arg
        elif opt == '--mp-is':
            mate_is = int(arg)
        elif opt == '--mp-dev':
            mate_is_dev = int(arg)
        elif opt == '--mp-orient':
            if arg == "fr" or arg == "rf":
                mate_orient = arg
            else:
                error("wrong value for --mp-orient option. Should be 'rf' or 'fr'")
        elif opt == '--min-k':
            min_k_rough = int(arg)
            k_2level_selection = False
        elif opt == '--max-k':
            max_k_rough = int(arg)
            k_2level_selection = False
        elif opt == '--step-k':
            step_k_rough = int(arg)
            k_2level_selection = False
        elif opt == '--mira-tmp-dir':
            mira_tmp_dir = arg
        elif opt == '--corrected':
            corrected = True
        elif opt == '--sc':
            single_cell = True
        elif opt == '--phred64':
            phred64 = True
        elif opt == '--assemblers':
            new_assemblers_to_run = arg.split(',')
            for name in new_assemblers_to_run:
                if not name in assemblers_to_run:
                    error("specified assembler name %s is not recognized! Use assemblers from this list: %s" % \
                          (name, str(assemblers_to_run)))
            assemblers_to_run = new_assemblers_to_run
        elif opt == '-t' or opt == "--threads":
            max_threads = int(arg)
        else:
            raise ValueError

    check_options()
    fill_in_common_params_subst_dict()

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    archive_dir = os.path.join(output_dir, ARCHIVE_DIR)
    best_dir = os.path.join(output_dir, BEST_DIR)
    if not os.path.isdir(archive_dir):
        os.makedirs(archive_dir)
    if not os.path.isdir(best_dir):
        os.makedirs(best_dir)

    log = logging.getLogger('runner')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    all_log_filename = os.path.join(archive_dir, "ALL_RUNNER" + LOG_SUFF)
    all_log_handler = logging.FileHandler(all_log_filename, mode='a')
    log.addHandler(all_log_handler)
    info("START RUNNER: " + " ".join(sys.argv) + "\n\n")

    # run dependent on K assemblers
    for (name, run_tool_function) in [('abyss', run_abyss), ('msrca', run_msrca), ('soap', run_soap),
        ('velvet', run_velvet), ('velvet-sc', run_velvet_sc), ('perga', run_perga), ('sga', run_sga), ('ray', run_ray)]:
        if not name in assemblers_to_run:
            continue
        if corrected and name == 'sga':
            continue
        cur_tool_log_filename = os.path.join(archive_dir, name + LOG_SUFF)
        cur_tool_log_handler = logging.FileHandler(cur_tool_log_filename, mode='a')
        log.addHandler(cur_tool_log_handler)
        iterate_K(name, run_tool_function)
        log.removeHandler(cur_tool_log_handler)

    # run independent on K assemblers
    for (name, run_tool_function) in [('spades', run_spades), ('idba', run_idba), ('sga', run_sga), ('cabog', run_cabog),
        ('mira', run_mira)]:
        if not name in assemblers_to_run:
            continue
        if not corrected and name == 'sga':
            continue
        cur_tool_log_filename = os.path.join(archive_dir, name + LOG_SUFF)
        cur_tool_log_handler = logging.FileHandler(cur_tool_log_filename, mode='a')
        log.addHandler(cur_tool_log_handler)
        run_K_independent_assembler(name, run_tool_function)
        log.removeHandler(cur_tool_log_handler)

    info("\n\n")
    info("FINISH RUNNER")


if __name__ == '__main__':
    main()
