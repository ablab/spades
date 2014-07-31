#!/usr/bin/python -O

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate coverage from raw file
import glob
import gzip
import sys
import os
import datetime
import getopt
from site import addsitedir
import logging
import shutil

from math import pow
from support import universal_sys_call, error
import options_storage

#profile = []
#insertions = {}
config = {}
#total_contigs
import os
import sys
import glob
import shutil
import support
import process_cfg
from site import addsitedir
from distutils import dir_util


def usage():
    sys.stderr.write('Mismatch Corrector. Simple post processing tool\n')
    sys.stderr.write('Usage: python' + str(sys.argv[0]) + '[options] -1 left_reads -2 right_reads -c contigs\n')
    sys.stderr.write('Or: python' + str(sys.argv[0]) + '[options] --12 mixed_reads -c contigs\n')
    sys.stderr.write('Or: python' + str(sys.argv[0]) + '[options] -s sam_file -c contigs\n')
    sys.stderr.write('Options:\n')
    sys.stderr.write('-t/--threads  <int>   threads number\n')
    sys.stderr.write('-o/--output-dir   <dir_name>  directory to store results\n')
    sys.stderr.write('-m/--mate-weight  <int>   weight for paired reads aligned properly. By default, equal to single reads weight (=1)\n')
    sys.stderr.write('--bwa <path>   path to bwa tool. Required if bwa is not in PATH\n')
    sys.stderr.write('--bowtie2    <path>  path to bowtie2 tool. Can be used instead of bwa. It is faster but provides a bit worse results\n')
    sys.stderr.write('  --use-quality use quality values as probabilities \n')
    sys.stderr.write('  --debug   save all intermediate files \n')
    sys.stderr.write('  --use-multiple-aligned  use paired reads with multiple alignment\n')
    sys.stderr.write('  --skip-masked   do not correct single \'N\' in contigs unless significant read support provided\n')
    sys.stderr.write('  --insert-size <int> estimation on insert size\n')
    sys.stderr.flush()


def parse_profile(args, log):
    global config

    long_options = "threads= sam-file= output-dir= bwa= contigs= mate-weight= dataset= split-dir= bowtie2= 12= insert-size= help debug use-quality use-multiple-aligned skip-masked".split()
    short_options = "1:2:o:s:S:c:t:m:q"

    reads1 = []
    reads2 = []
    reads_mixed = []

    def check_file(f, type, log):
        f = os.path.abspath(f)
        if not os.path.isfile(f):
            error("file with %s (%s) doesn't exist! " % (type, f), log)
        return f

    options, contigs_fpaths = getopt.gnu_getopt(args, short_options, long_options)
    for opt, arg in options:
        if opt in ('-o', "--output-dir"):
            config["output_dirpath"] = os.path.abspath(arg)
            config["make_latest_symlink"] = False
        if opt in ('-c', "--contigs"):
            config["contigs"] = check_file(arg, "contigs", log)
        if opt == '-1':
            reads1.append(check_file(arg, 'left reads', log))
        if opt == '-2':
            reads2.append(check_file(arg, 'right reads', log))
#            config["reads2"] = os.path.abspath(arg)
#            if not os.path.exists(config["reads2"]):
#                log.info("FILE WITH READS DOES NOT EXIST!")
#                usage()
#                sys.exit(1)
        if opt == '--12':
            reads_mixed.append(check_file(arg, 'interleaved left and right reads', log))
#            config["reads_mixed"] = os.path.abspath(arg)
#            if not os.path.exists(config["reads_mixed"]):
#                log.info("FILE WITH READS DOES NOT EXIST!")
#                usage()
#                sys.exit(1)
        if opt == "--bwa":
            config["bwa"] = os.path.abspath(arg)
        if opt in ('-t', "--threads"):
            config["t"] = int(arg)
        if opt in ('-m', "--mate-weight"):
            config["mate_weight"] = float(arg)
        if opt in ('-s', "--sam-file"):
            config["sam_file"] = os.path.abspath(arg)
        if opt == "--dataset":
            config["dataset"] = os.path.abspath(arg)
        if opt in ('-S', "--split-dir"):
            config["split_dir"] = os.path.abspath(arg)
        if opt in ('-q', "--use-quality"):
            config["use_quality"]= 1
        if opt == "--bowtie2":
            if arg != "bowtie2":
                arg = os.path.abspath(arg)
            config["bowtie2"]= arg
        if opt == "--debug":
            config["debug"] = 1
        if opt in ("--use-multiple-aligned"):
            config["use_multiple_aligned"] = 1
        if opt == "--skip-masked":
            config["skip_masked"] = 1
        if opt == "--insert-size":
            config["insert_size"] = int(arg)

    if len(reads1) != len(reads2):
        error("the number of files with left paired reads is not equal to the"
              " number of files with right paired reads!", log)

    if reads1:
        config["reads1"] = reads1
    if reads2:
        config["reads2"] = reads2
    if reads_mixed:
        config["reads_mixed"] = reads_mixed

    work_dir = os.path.join(config["output_dirpath"], "tmp")
    config["work_dir"] = work_dir
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
#    os.system ("mkdir -p " + work_dir)
#    os.system("rm -rf " + work_dir + "/*")


def init_config():
    global config
    config = {}

    now = datetime.datetime.now()
    config["output_dirpath"] = "corrector.output." + now.strftime("%Y.%m.%d_%H.%M.%S")+"/"
    config["bwa"] = "bwa"
    config["t"] = int(4)
    config["mate_weight"] = float(1)
    config["use_quality"] = 0
    config["debug"] = 0
    config["use_multiple_aligned"] = 0
    config["skip_masked"] = 0
    config["insert_size"] = int(400)


def process_contig(files):
    log = logging.getLogger('spades')
    samfilename = files[0]
    contig_file = files[1]
    #contigs are always .fasta
    contig_name = contig_file[0:-6]
    logFileName =   contig_name + '.stdout'
    logFile = open(logFileName, 'w')
    ntime = datetime.datetime.now()
    starttime = ntime
  #  os.system ('./build/release/bin/corrector ' + samfilename + ' ' + contig_file)

    logFile.write(ntime.strftime("%Y.%m.%d_%H.%M.%S") + ": All done. ")
    stime = ntime - starttime
    logFile.write("Time spent: " + str(stime) + "\n")

    logFile.write("Finished processing "+ str(contig_file)) # + ". Used " + str(total_reads) + " reads.\n")
    #logFile.write("replaced: " + str(replaced) + " deleted: "+ str(deleted) +" inserted: " + str(inserted) +'\n')
    logFile.close()

#    return inserted, replaced

def prepare_config_corr(filename, cfg, ext_python_modules_home):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml
    #print "dumping contigs to " + filename
    data = pyyaml.load(open(filename, 'r'))
    data["dataset"] = cfg.dataset_yaml_filename
    data["output_dir"] = cfg.output_dir
    data["work_dir"] = cfg.output_dir + '/tmp'
    #data["hard_memory_limit"] = cfg.max_memory
    data["max_nthreads"] = cfg.max_threads
    data["bwa"] = cfg.bwa
    file_c = open(filename, 'w')
    pyyaml.dump(data, file_c)
    file_c.close()


def main(args, joblib_path, log=None, config_file=None):

    if len(args) < 1:
        usage()
        sys.exit(0)
    addsitedir(joblib_path)

    init_config()
    parse_profile(args, log)

    if not log:
        log = logging.getLogger('spades')
        log.setLevel(logging.DEBUG)

        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.DEBUG)
        log.addHandler(console)

        log_filename = os.path.join(config["output_dirpath"], "corrector.log")
        log_handler = logging.FileHandler(log_filename, mode='w')
        log.addHandler(log_handler)

    log.info("Config: " + str(config))
    if "split_dir" not in config:
#        print "no split dir, looking for sam file"
        '''        if "sam_file" not in config:
            log.info("no sam file, running aligner")
            run_aligner(log)
        else:
            log.info("sam file was found")
            tmp_sam_file_path = os.path.join(config["work_dir"], "tmp.sam")
            shutil.copy2(config["sam_file"], tmp_sam_file_path) # Note: shutil.copy2 is similar to the Unix command cp -p
            #os.system("cp -p "+ config["sam_file"] +" " + config["work_dir"]+"tmp.sam")
            config["sam_file"] = tmp_sam_file_path '''
        path_to_bin = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../bin/corrector')
        path_to_config = os.path.join(os.path.dirname(os.path.realpath(__file__)) , '../../configs/corrector/corrector.info.template')
        if config_file:
            path_to_config = config_file
       # config["output_dirpath"] += "/mismatch_corrector_tmp"
#        print config["output_dirpath"] + " output_dirpath"
#        print path_to_config
        run_str = path_to_bin + ' ' + path_to_config + ' ' + config["contigs"]
#        print run_str
        os.system (run_str)
    #    now = datetime.datetime.now()
    #    res_directory = "corrector.output." + now.strftime("%Y.%m.%d_%H.%M.%S")+"/"



def run_corrector(corrected_dataset_yaml_filename, configs_dir, execution_home, cfg,
                ext_python_modules_home, log, args):
    addsitedir(ext_python_modules_home)
    if sys.version.startswith('2.'):
        import pyyaml2 as pyyaml
    elif sys.version.startswith('3.'):
        import pyyaml3 as pyyaml
    cfg.dataset_yaml_filename = corrected_dataset_yaml_filename
    #print cfg.dataset_yaml_filename + " yaml with all_libs"
    #print configs_dir + " configs dir"
    #print execution_home + " execution home"
    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    dir_util.copy_tree(os.path.join(configs_dir, "corrector"), dst_configs, preserve_times=False)
    cfg_file_name = os.path.join(dst_configs, "corrector.info")
    # removing template configs
    for root, dirs, files in os.walk(dst_configs):
        for cfg_file in files:
            cfg_file = os.path.join(root, cfg_file)
            if cfg_file.endswith('.template'):
                if os.path.isfile(cfg_file.split('.template')[0]):
                    os.remove(cfg_file)
                else:
                    os.rename(cfg_file, cfg_file.split('.template')[0])

    cfg.tmp_dir = support.get_tmp_dir(prefix="corrector_")

    prepare_config_corr(cfg_file_name, cfg, ext_python_modules_home)
    binary_name = "corrector"

    command = [os.path.join(execution_home, binary_name),
               os.path.abspath(cfg_file_name)]

    log.info("\n== Running contig polishing tool: " + ' '.join(command) + "\n")
 #   support.sys_call(command, log)
 #   if not os.path.isfile(corrected_dataset_yaml_filename):
 #       support.error("read error correction finished abnormally: " + corrected_dataset_yaml_filename + " not found!")
#    corrected_dataset_data = pyyaml.load(open(corrected_dataset_yaml_filename, 'r'))

    is_changed = False

    log.info("\n== Dataset description file was created: " + cfg_file_name + "\n")
    main(args, ext_python_modules_home, log, cfg_file_name)



if __name__ == '__main__':
    joblib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../ext/src/python_libs')
    main(sys.argv[1:], joblib_path)
    


