#!/usr/bin/python

import os
import sys
import subprocess
import itertools
import matplotlib
import shutil

########################################################################	
### CONFIG & CHECKS 
########################################################################

BUILD_PATH = os.path.join(os.path.abspath(sys.path[0]), '../../../build/ext/allpaths')
WORK_DIR = os.path.join(BUILD_PATH, 'temp')
DATA_DIR = os.path.join(WORK_DIR, 'genome/data')

########################################################################
def parse_config_line(config_file, line, value_name, value):
    if line.startswith(value_name) and len(line.split()) > 1:
        value=line.split()[1]
        line = ""
        while line.strip() == "":
            line = config_file.readline()
            if not line:
                print "incorrect config. Exitting.."
                sys.exit(0)
    return line, value    

########################################################################

if len(sys.argv) != 2:
    print 'ALLPATHS runner'
    print 'Usage: ', sys.argv[0], ' CONFIG_FILE (See README)'
    sys.exit(0)

# check config
if not os.path.isfile(sys.argv[0]):
    print "config file not found. Exitting.."
    sys.exit(0)

config_file = open(sys.argv[1])

frag_lib=""
frag_is=0
frag_dev=0

jump_lib=""
jump_is=0
jump_dev=0

contigs_for_sim=""
ploidy=1
output_contigs=""

line = config_file.readline()
line, frag_lib  = parse_config_line(config_file, line, "frag_lib", frag_lib)
line, frag_is   = parse_config_line(config_file, line, "insert_size", frag_is)
frag_is=int(frag_is)
line, frag_dev  = parse_config_line(config_file, line, "stddev", frag_dev)
frag_dev=int(frag_dev)

line, jump_lib  = parse_config_line(config_file, line, "jump_lib", jump_lib)
line, jump_is   = parse_config_line(config_file, line, "insert_size", jump_is)
jump_is=int(jump_is)
line, jump_dev  = parse_config_line(config_file, line, "stddev", jump_dev)
jump_dev=int(jump_dev)

line, contigs_for_sim = parse_config_line(config_file, line, "contigs_for_sim", contigs_for_sim)
line, ploidy          = parse_config_line(config_file, line, "ploidy", ploidy)           
ploidy=int(ploidy)

if line.startswith("output_contigs"):
    output_contigs=line.split()[1]

frag_lib        = os.path.expanduser(frag_lib)
jump_lib        = os.path.expanduser(jump_lib)
output_contigs  = os.path.expanduser(output_contigs)
contigs_for_sim = os.path.expanduser(contigs_for_sim)

if (output_contigs == "" or (frag_lib == "" and jump_lib == "") or ( (frag_lib == "" or jump_lib == "") and contigs_for_sim == "" ) ):
    print "incorrect config. Exitting.."
    sys.exit(0)

##############################
# copy sources and make
subprocess.call(['sh', 'run_allpaths_helper.sh', 'make']) 

# check if simulation is needed
jump_inward = False
if (frag_lib == "" or jump_lib == ""):
    if not os.path.isfile(contigs_for_sim):
        print "contigs for simulation not found. Exitting.."
        sys.exit(0)

    wgsim_path = os.path.join(BUILD_PATH, 'wgsim_src/wgsim')  
    sim_lib_dir = os.path.dirname(output_contigs)
    insert_size = frag_is
    stddev      = frag_dev
    if (jump_lib == ""):
        insert_size = jump_is
        stddev      = jump_dev          

    print "== simulating reads (" + sim_lib_dir + "/out.read1.fq and out.read2.fq) =="
    sys.stdout.flush()
    sys.stderr.flush()
    logfile_out = open(os.path.dirname(wgsim_path) + '/log.txt', 'w')
    logfile_err = open(os.path.dirname(wgsim_path) + '/log.err', 'w')
    subprocess.call([wgsim_path, contigs_for_sim, sim_lib_dir + '/out.read1.fq', sim_lib_dir + '/out.read2.fq', '-d', str(insert_size), '-s', str(stddev), 
        '-N', '1000000', '-1', '100', '-2', '100', '-e', '0', '-r', '0', '-R', '0', '-X', '0'], stdout=logfile_out, stderr=logfile_err)
    logfile_out.close()
    logfile_err.close()
    print "== simulating finished =="

    if (frag_lib == ""):
        frag_lib = sim_lib_dir + "/out.read?.fq"       
    else:
        jump_lib = sim_lib_dir + "/out.read?.fq"    
        jump_inward = True   

##############################
# preparing data for allpaths
if not os.path.isdir(WORK_DIR):
    os.mkdir(WORK_DIR)

in_groups = open(WORK_DIR + '/in_groups.csv', 'w')
in_groups.write("library_name, group_name, file_name\n")
in_groups.write("my_frag,      frags,      " + frag_lib + "\n")
in_groups.write("my_jump,      jumps,      " + jump_lib)
in_groups.close()

in_libs = open(WORK_DIR + '/in_libs.csv', 'w')
in_libs.write("library_name, project_name, organism_name,     type, paired, frag_size, frag_stddev, insert_size, insert_stddev, read_orientation, genomic_start, genomic_end\n")
in_libs.write("my_frag,           my_proj,      organism, fragment,      1, " + str(frag_is) + ", "  + str(frag_dev) + ",  ,  ,           inward,             0,           0\n")
in_libs.write("my_jump,           my_proj,      organism,  jumping,      1,  ,  , " + str(jump_is) + ", "  + str(jump_dev) + ", ")
if jump_inward:
    in_libs.write("inward")
else:
    in_libs.write("outward")
in_libs.write(",             0,           0")
in_libs.close()

prepare = open(WORK_DIR + '/prepare.sh', 'w')
prepare.write("#!/bin/sh\n")
prepare.write("# ALLPATHS-LG needs 100 MB of stack space.  In 'csh' run 'limit stacksize 100000'.\n")
prepare.write("ulimit -s 100000\n")
prepare.write("\n")
prepare.write("mkdir -p " + DATA_DIR + "\n")
prepare.write("\n")
prepare.write("PrepareAllPathsInputs.pl\\\n")
prepare.write(" DATA_DIR=" + DATA_DIR + "\\\n")
prepare.write(" PLOIDY=" + str(ploidy) + "\\\n")
prepare.write(" IN_GROUPS_CSV=" + WORK_DIR + '/in_groups.csv' + "\\\n")
prepare.write(" IN_LIBS_CSV=" + WORK_DIR + '/in_libs.csv' + "\\\n")
prepare.write(" OVERWRITE=True\\\n")
prepare.write(" PHRED_64=0\n")
prepare.close()

subprocess.call(['sh', 'run_allpaths_helper.sh', 'prepare'])

##############################
# run assembler

assemble = open(WORK_DIR + '/assemble.sh', 'w')
assemble.write("#!/bin/sh\n")
assemble.write("# ALLPATHS-LG needs 100 MB of stack space.  In 'csh' run 'limit stacksize 100000'.\n")
assemble.write("ulimit -s 100000\n")
assemble.write("\n")
assemble.write("RunAllPathsLG \\\n")
assemble.write(" PRE=" + WORK_DIR + "\\\n")
assemble.write(" REFERENCE_NAME=genome\\\n")
assemble.write(" DATA_SUBDIR=data\\\n")
assemble.write(" RUN=run\\\n")
assemble.write(" SUBDIR=test\\\n")
assemble.write(" TARGETS=standard\\\n")
assemble.write(" OVERWRITE=True\n")
assemble.close()

subprocess.call(['sh', 'run_allpaths_helper.sh', 'assemble'])

##############################
# moving contigs to final place
contigs = DATA_DIR + '/run/test/final.contigs.fasta'
scaffold = DATA_DIR + '/run/test/final.assembly.fasta'
if os.path.isfile(contigs):
    shutil.move(contigs, output_contigs)
else:
    print "FINAL CONTIGS NOT FOUND! See logs for errors"
if os.path.isfile(scaffold):
    basename, extension = os.path.splitext(output_contigs)
    output_scaffold = basename + '_scaffold' + extension
    shutil.move(scaffold, output_scaffold)
else:
    print "FINAL SCAFFOLD NOT FOUND! See logs for errors"

##############################
# clearing temp data
if os.path.isdir(DATA_DIR):
    shutil.rmtree(DATA_DIR)

