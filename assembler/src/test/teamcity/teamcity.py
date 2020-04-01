#!/usr/bin/python

############################################################################
# Copyright (c) 2015-2017 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

#script for testing SPAdes
#provide a path to .info file

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

sys.path.append('./src/test/teamcity/')
import teamcity_support
from teamcity_support import *

### main ###
try:
    if len(sys.argv) == 1:
        command = 'python {} -h'.format(sys.argv[0])
        subprocess.call(command, shell=True)
        sys.exit(1)

    sys.stderr = sys.stdout
    exit_code = 0
    args = parse_args()
    dataset_info = load_info(args.info)
    working_dir = os.getcwd()
    output_dir = create_output_dir(args, dataset_info)
    save_run_info(args, output_dir)

    #compile
    log.log("##teamcity[blockOpened name='build' description='Building SPAdes']")
    ecode = compile_spades(args, dataset_info, working_dir)
    if ecode != 0:
        log.err("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)
    log.log("##teamcity[blockClosed name='build']")

    #run spades
    log.log("##teamcity[blockOpened name='run' description='Running SPAdes']")
    spades_dir = working_dir
    if args.spades_path:
        spades_dir = args.spades_path
        log.log("Different spades.py path specified: " + spades_dir)
    spades_cmd = make_spades_cmd(args, dataset_info, spades_dir, output_dir)
    #log.log("Launching: " + spades_cmd)

    ecode = os.system(spades_cmd)
    if ecode != 0:
        log.err("SPAdes finished abnormally with exit code " + str(ecode))
        sys.exit(4)
    log.log("##teamcity[blockClosed name='run']")

    #reads quality
    if 'reads_quality_params' in dataset_info.__dict__:
        log.log("##teamcity[blockOpened name='reads quality' description='Asessing reads quality']")
        exit_code = run_reads_assessment(dataset_info, working_dir, output_dir)
        log.log("##teamcity[blockClosed name='reads quality']")


    #QUAST
    rewrite_latest = True
    contigs = get_contigs_list(args, dataset_info)
    if 'quast_params' in dataset_info.__dict__:
        log.log("##teamcity[blockOpened name='quast' description='Running QUAST']")
        ecode = quast_analysis(contigs, dataset_info, output_dir)
        if ecode != 0:
            rewrite_latest = False
            log.err("QUAST analysis did not pass, exit code " + str(ecode))
            exit_code = ecode
        log.log("##teamcity[blockClosed name='quast']")


    #etalon saves
    if 'etalon_saves' in dataset_info.__dict__:
        log.log("Comparing etalon saves now")
        log.log("##teamcity[blockOpened name='saves' description='Comparing etalon saves']")
        ecode = os.system(os.path.join(spades_dir, "./src/test/teamcity/detect_diffs.sh") + " " + output_dir + " " + dataset_info.etalon_saves)
        if ecode != 0:
            rewrite_latest = False
            log.err("Comparing etalon saves did not pass, exit code " + str(ecode))
            exit_code = 12
        log.log("##teamcity[blockClosed name='saves']")

    #compare misassemblies
    log.log("##teamcity[blockOpened name='misassemblies' description='Comparing misassemblies']")
    contig_storage_dir = get_contigs_storage_dir(args, dataset_info)
    ecode, rewrite = compare_misassemblies(contigs, dataset_info, contig_storage_dir, output_dir)
    rewrite_latest = rewrite_latest and rewrite
    if ecode != 0:
        log.err('Failed to compare misassemblies')
    log.log("##teamcity[blockClosed name='misassemblies']")

    #save contigs to storage
    contigs = get_contigs_list(args, dataset_info, True)
    save_contigs(args, output_dir, contig_storage_dir, contigs, rewrite_latest)

    #save quast report as build artifact
    artifact_dir = os.path.join(working_dir, "quast_reports")
    save_quast_report(contigs, dataset_info, contig_storage_dir, output_dir, artifact_dir)

    sys.exit(exit_code)

except SystemExit:
    raise

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(239)
