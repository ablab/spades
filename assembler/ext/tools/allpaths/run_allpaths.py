#!/usr/bin/python

import os
import sys
import subprocess
import shutil
import glob

########################################################################	
### CONFIG & CHECKS 
########################################################################

BUILD_PATH = os.path.join(os.path.abspath(sys.path[0]), '../../../build/ext/allpaths')
WORK_DIR = os.path.join(BUILD_PATH, 'temp')
DATA_DIR = os.path.join(WORK_DIR, 'genome/data')

########################################################################

def parse_config_line(config_file, line, value_name, value):
    if line.startswith(value_name) and len(line.split(value_name)) > 1:
        value = line.split(value_name)[1].split(";")[0].strip()                
        line = ""
        while line.strip() == "":
            line = config_file.readline()
            if not line:
                print "incorrect config. Exitting.."
                sys.exit(0)
    return line, value  

def get_quality_encoding(filename):
    fastq_file = ""
    if filename.endswith(".gz"):
        import gzip
        fastq_file = gzip.open(filename, 'r')    
    else:
        fastq_file = open(filename, 'r')
    i = 0
    for line in fastq_file:
        if i == 3:
            for c in line:
                if ord(c) < 59:
                    fastq_file.close()
                    return "phred+33"
                if ord(c) > 74:
                    fastq_file.close()
                    return "phred+64"
        i = (i + 1) % 4
    fastq_file.close()

    return "not_defined"

def convert_sim_lib(sim_lib_read1, sim_lib_read2):
    print "== converting simulated library to Phred+64 =="
    import fileinput
    for lines in fileinput.input(sim_lib_read1, inplace = 1): ## edit file in place
        lines = lines.replace("I","h")
        print lines.strip()
    for lines in fileinput.input(sim_lib_read2, inplace = 1): ## edit file in place
        lines = lines.replace("I","h")
        print lines.strip()
    print "== converting finished =="
    return

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

if ( frag_lib != "" and len(glob.glob(frag_lib)) == 0 ):
    print "fragment library (" + frag_lib + ") not found. Exitting.."
    sys.exit(0)
if ( jump_lib != "" and len(glob.glob(jump_lib)) == 0 ):
    print "jumping library (" + jump_lib + ") not found. Exitting.."
    sys.exit(0)
if ( contigs_for_sim != "" and not os.path.isfile(contigs_for_sim) ):
    print "contigs for simulation (" + contigs_for_sim + ") not found. Exitting.."
    sys.exit(0)

if (output_contigs == "" or (frag_lib == "" and jump_lib == "") or ( (frag_lib == "" or jump_lib == "") and contigs_for_sim == "" ) ):
    print "incorrect config. Exitting.."
    sys.exit(0)

##############################
# copy sources and make
subprocess.call(['sh', 'helper_run_allpaths.sh', 'make']) 

# check if simulation is needed
frag_simulated = False
jump_simulated = False
sim_lib_read1  = ""
sim_lib_read2  = ""
if (frag_lib == "" or jump_lib == ""):
    if not os.path.isfile(contigs_for_sim):
        print "contigs for simulation not found. Exitting.."
        sys.exit(0)

    wgsim_path = os.path.join(BUILD_PATH, 'wgsim_src/wgsim')  
    sim_lib_dir = os.path.dirname(output_contigs)
    if not os.path.isdir(sim_lib_dir):
        os.makedirs(sim_lib_dir)

    insert_size = frag_is
    stddev      = frag_dev
    if (jump_lib == ""):
        insert_size = jump_is
        stddev      = jump_dev          

    sim_lib_read1 = sim_lib_dir + "/out.read1.fq"
    sim_lib_read2 = sim_lib_dir + "/out.read2.fq"
    print "== simulating reads (" + sim_lib_read1 + " and " + sim_lib_read2 + ") =="
    sys.stdout.flush()
    sys.stderr.flush()
    logfile_out = open(os.path.dirname(wgsim_path) + '/sim.log', 'w')
    logfile_err = open(os.path.dirname(wgsim_path) + '/sim.err', 'w')
    subprocess.call([wgsim_path, contigs_for_sim, sim_lib_read1, sim_lib_read2, '-d', str(insert_size), '-s', str(stddev), 
        '-N', '1000000', '-1', '100', '-2', '100', '-e', '0', '-r', '0', '-R', '0', '-X', '0'], stdout=logfile_out, stderr=logfile_err)
    logfile_out.close()
    logfile_err.close()
    print "== simulating finished =="

    if (frag_lib == ""):
        frag_lib = sim_lib_dir + "/out.read?.fq"  
        frag_simulated = True
    else:
        jump_lib = sim_lib_dir + "/out.read?.fq"    
        jump_simulated = True   

##############################
# preparing data for allpaths
if not os.path.isdir(WORK_DIR):
    os.makedirs(WORK_DIR)

in_groups = open(WORK_DIR + '/in_groups.csv', 'w')
in_groups.write("library_name, group_name, file_name\n")
in_groups.write("my_frag,      frags,      " + frag_lib + "\n")
in_groups.write("my_jump,      jumps,      " + jump_lib)
in_groups.close()

in_libs = open(WORK_DIR + '/in_libs.csv', 'w')
in_libs.write("library_name, project_name, organism_name,     type, paired, frag_size, frag_stddev, insert_size, insert_stddev, read_orientation, genomic_start, genomic_end\n")
in_libs.write("my_frag,           my_proj,      organism, fragment,      1, " + str(frag_is) + ", "  + str(frag_dev) + ",  ,  ,           inward,             0,           0\n")
in_libs.write("my_jump,           my_proj,      organism,  jumping,      1,  ,  , " + str(jump_is) + ", "  + str(jump_dev) + ", ")
if jump_simulated:
    in_libs.write("inward")
else:
    in_libs.write("outward")
in_libs.write(",             0,           0")
in_libs.close()

##############################
# checking quality encoding
frag_quality = "phred+33"
jump_quality = "phred+33"
if not frag_simulated:
    frag_1 = glob.glob(frag_lib)[0]
    if not os.path.isfile(frag_1):
        print "fragment library (" + frag_1 + ") not found. Exitting.."
        sys.exit(0)
    frag_quality = get_quality_encoding(frag_1)
if not jump_simulated:
    jump_1 = glob.glob(jump_lib)[0]
    if not os.path.isfile(jump_1):
        print "jumping library (" + jump_1 + ") not found. Exitting.."
        sys.exit(0)
    jump_quality = get_quality_encoding(jump_1)
    
if jump_quality == "not_defined" or frag_quality == "not_defined":
    print "can't define quality score of one or more libraries. Exitting.."
    sys.exit(0)

PHRED_64 = 0
if jump_quality != frag_quality:
    if not jump_simulated and not frag_simulated:
        print "different quality score in jumping and fragment library. Exitting.."
        sys.exit(0)
    else: # simulated library always in Phred+33 format, so we should convert it to Phred+64
        convert_sim_lib(sim_lib_read1, sim_lib_read2)
        PHRED_64 = 1
else:
    if frag_quality == "phred+33":
        PHRED_64 = 0
    else:   
        PHRED_64 = 1           

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
prepare.write(" PHRED_64=" + str(PHRED_64) + "\n")
prepare.close()

subprocess.call(['sh', 'helper_run_allpaths.sh', 'prepare'])

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

subprocess.call(['sh', 'helper_run_allpaths.sh', 'assemble'])

##############################
# moving contigs to final place
contigs = DATA_DIR + '/run/ASSEMBLIES/test/final.contigs.fasta'
scaffold = DATA_DIR + '/run/ASSEMBLIES/test/final.assembly.fasta'
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

