#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import shutil
import re
import getopt
import datetime
import subprocess

###################################################################

sys.path.append(os.path.join(os.path.abspath(sys.path[0]), 'conversion'))
sys.path.append(os.path.join(os.path.abspath(sys.path[0]), 'stat'))

bowtie_path  = os.path.join(os.path.abspath(sys.path[0]), '../../../../external_tools/bowtie-0.12.7')
bowtie_build = os.path.join(bowtie_path, "bowtie-build")
bowtie       = os.path.join(bowtie_path, "bowtie")

tmp_folder = "tmp"
output_dir = "results_" + datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
thread_num = 16
bin_size = 1
kmer = 1
make_latest_symlink = True
reference = ""
max_is = 1000000000
skip_trimming = False

###################################################################

long_options = "output-dir= reference= thread-num= bin-size= kmer-size= max-is= skip-trimming".split()
short_options = "o:r:t:b:k:x:s"

def usage():
    print 'Estimation reads quality'
    print 'Usage:', sys.argv[0], ' [options described below] datasets description-file(s)'
    print ""
    print "Options with parameters:"
    print "-r\t--reference\tFile with reference genome (Mandatory parameter)"
    print "-o\t--output-dir\tDirectory to store all result files"
    print "-t\t--thread-num\tMax number of threads (default is " + str(thread_num) + ")"
    print "-x\t--max-is\tMaximal inser size (default is none)"
    print "-s\t--skip-trimming\tSkip N-trimming for speed-up"
    
def check_file(f):
    if not os.path.isfile(f):
        print "Error - file not found:", f
        sys.exit(2)
    return f

#####
print "\n", sys.argv[0], "is DEPRECATED!\nUse reads_quality.py --paired-mode instead!", "\n"
sys.exit(1)
#####

try:
    options, datasets = getopt.gnu_getopt(sys.argv[1:], short_options, long_options)
except getopt.GetoptError, err:
    print str(err)
    print ""
    usage()
    sys.exit(1)

for opt, arg in options:
    if opt in ('-o', "--output-dir"):
        output_dir = arg
        make_latest_symlink = False  
    elif opt in ('-r', "--reference"):
        reference = arg
    elif opt in ('-t', "--thread-num"):
        thread_num = int(arg)
        if thread_num < 1:
            thread_num = 1 
    elif opt in ('-x', "--max-is"):
        max_is = int(arg)
        if max_is < 0:
            max_is = 1000000000
    elif opt in ('-s', "--skip-trimming"):
        skip_trimming = True
    else:
        raise ValueError

for d in datasets:
    check_file(d)    

if not datasets:
    print "no datasets"
    usage()
    sys.exit(1)   

if not reference:
    print 'no ref'
    usage()
    sys.exit(1)   

###################################################################

def get_full_path(dataset, rel_path):
    return os.path.abspath(os.path.join(os.path.dirname(dataset), rel_path))

def ungzip_if_needed(filename, output_folder, force_copy = False):
    file_basename, file_extension = os.path.splitext(filename)
    if file_extension == ".gz":
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        ungzipped_filename = os.path.join(output_folder, os.path.basename(file_basename))
        ungzipped_file = open(ungzipped_filename, 'w')
        subprocess.call(['gunzip', filename, '-c'], stdout=ungzipped_file)
        ungzipped_file.close()
        filename = ungzipped_filename    
    elif force_copy:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        shutil.copy(filename, output_folder)
        filename = os.path.join(output_folder, os.path.basename(filename))
    return filename

###################################################################

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if make_latest_symlink:
    latest_symlink = 'latest'
    if os.path.islink(latest_symlink):
        os.remove(latest_symlink)
    os.symlink(output_dir, latest_symlink)

datasets_dict = dict()

print("Analyzing datasets")
for dataset in datasets:

    basename = os.path.splitext(os.path.basename(dataset))[0]
    cur_key = basename
    i = 1
    while datasets_dict.has_key(cur_key):
        cur_key = basename + "_" + str(i)

    cur_reads = []
    for line in open(dataset, 'r'):
        if line.startswith("paired_reads") or line.startswith("single_reads"):
            line = line.replace('"', '')
            reads = line.split()[1:]
            for read in reads:
                cur_reads.append(get_full_path(dataset, read))            
    
    if len(cur_reads) == 0:
        print("  " + dataset + " was skipped because it contains no reads")    
        continue
        
    datasets_dict[cur_key] = cur_reads   
    print("  " + dataset + " ==> " + cur_key)

if len(datasets_dict.keys()) == 0:
    print("Can't continue estimation - all datasets were skipped")
    sys.exit(1)

###################################################################

report_dict = {"header" : ["Dataset"]}
for dataset in datasets_dict.iterkeys():
    report_dict[dataset] = [dataset]

tmp_folder = os.path.join(output_dir, tmp_folder)
if not os.path.exists(tmp_folder):
    os.makedirs(tmp_folder)

if not skip_trimming:
    print("Unpacking data (if needed) to temporary folder (" + tmp_folder + ") and N-trimming")
else:
    print("Unpacking data (if needed) to temporary folder (" + tmp_folder + ")")

print("  reference...")
reference = ungzip_if_needed(reference, tmp_folder)

# TODO fastA analysis (we should convert all in fasta if there is at least one file in fasta)
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    ungzipped_reads = []
    for read in datasets_dict[dataset]:    
        copied_read = ungzip_if_needed(read, os.path.join(tmp_folder, dataset), True)
        if not skip_trimming:
            import trim_ns
            trim_ns.trim_file(copied_read, copied_read)
        ungzipped_reads.append(copied_read)
    datasets_dict[dataset] = ungzipped_reads

# creating index
index_folder = os.path.join(tmp_folder, "index")
if not os.path.exists(index_folder):
    os.makedirs(index_folder)
index_name   = os.path.splitext(os.path.basename(reference))[0]
index        = os.path.join(index_folder, index_name)
print("Creating index " + index)
index_log = open(os.path.join(output_dir, "index.log"),'w')
index_err = open(os.path.join(output_dir, "index.err"),'w')
subprocess.call([bowtie_build, reference, index], stdout=index_log, stderr=index_err)
index_log.close()
index_err.close()

# bowtie-ing
print("Aligning (ignoring reads with multiple possible aligment)")
report_dict["header"] += ["Total reads", "Uniquely aligned reads", "Unaligned reads", "Non-niquely aligned reads"]
total_reads = {}
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = open(os.path.join(output_dir, dataset + ".log"),'w')
    align_err = open(os.path.join(output_dir, dataset + ".err"),'w') 
    reads_string = reduce(lambda x, y: x + ',' + y, datasets_dict[dataset])   
    subprocess.call([bowtie, '-c', '-q', '-m', '1', '--suppress', '6,7,8', index, '-p', str(thread_num), reads_string], stdout=align_log, stderr=align_err)
    align_log.close()
    align_err.close() 

    align_err = open(os.path.join(output_dir, dataset + ".err"),'r') 
    suppressed_added = False
    for line in align_err:
        if line.startswith("# reads processed") or line.startswith("# reads with at least one") or line.startswith("# reads that failed"):
            report_dict[dataset].append( (line.split(':')[1]).strip() )
        elif line.startswith("# reads with alignments suppressed due to"):
            report_dict[dataset].append( (line.split(':')[1]).strip() )
            suppressed_added = True
        if line.startswith("# reads processed"):
            total_reads[dataset] = int((line.split(':')[1]).strip())

    align_err.close()       

    if not suppressed_added:
        report_dict[dataset].append( "0 (0.00%)" )

# raw-single    
print("Parsing Bowtie log")
import raw_single
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = os.path.join(output_dir, dataset + ".log")
    raw_file  = os.path.join(output_dir, dataset + ".raw")
    raw_single.raw_single(align_log, raw_file)

# get length of reference
ref_len = 0
for line in open(reference):
    if line[0] != '>':       
        ref_len += len(line.strip())

# coverage # python reads_utils/stat/coverage.py ec.raw ec.cov 4639675 1000
print("Analyzing coverage")
report_dict["header"] += ["Genome mapped (%)"]
gaps_dict = {}
import coverage
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    raw_file  = os.path.join(output_dir, dataset + ".raw")
    cov_file  = os.path.join(output_dir, dataset + ".cov")
    cov = coverage.coverage(raw_file, cov_file, ref_len, 1, 1)
    gaps_file  = os.path.join(output_dir, dataset + ".gaps")
    chunks_file  = os.path.join(output_dir, os.path.splitext(os.path.basename(reference))[0] + "gaps_" + dataset + ".fasta")
    gaps_dict[dataset] = coverage.analyze_gaps(cov_file, gaps_file, reference, chunks_file, kmer)

    report_dict[dataset].append( str(cov * 100) )

# is form logs    
print("Retaining insert size")
report_dict["header"] += ["Read length", "FR read pairs", "Insert size (deviation)", "RF read pairs", "Insert size (deviation)", "FF read pairs", "Insert size (deviation)", "One uniquely aligned read in pair", "Suppressed due to insert size limit"]
import is_from_single_log
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = os.path.join(output_dir, dataset + ".log")
    stat = is_from_single_log.stat_from_log(align_log, max_is)
    stat[1]["FR"].write_hist(os.path.join(output_dir, dataset + "_FR.is"))
    stat[1]["RF"].write_hist(os.path.join(output_dir, dataset + "_RF.is"))
    stat[1]["FF"].write_hist(os.path.join(output_dir, dataset + "_FF.is"))

    read_pairs = total_reads[dataset] / 2

    report_dict[dataset].append( str(stat[0]) )

    report_dict[dataset].append( str(stat[1]["FR"].count) + " (" + str(round( 100.0 * float(stat[1]["FR"].count) / float(read_pairs), 2) ) + "%)" )
    report_dict[dataset].append( str(round(stat[1]["FR"].mean, 2)) + " (" + str(round(stat[1]["FR"].dev, 2)) + ")"  )

    report_dict[dataset].append( str(stat[1]["RF"].count) + " (" + str(round( 100.0 * float(stat[1]["RF"].count) / float(read_pairs), 2) ) + "%)" )
    report_dict[dataset].append( str(round(stat[1]["RF"].mean, 2)) + " (" + str(round(stat[1]["RF"].dev, 2)) + ")"  )

    report_dict[dataset].append( str(stat[1]["FF"].count) + " (" + str(round( 100.0 * float(stat[1]["FF"].count) / float(read_pairs), 2) ) + "%)" )
    report_dict[dataset].append( str(round(stat[1]["FF"].mean, 2)) + " (" + str(round(stat[1]["FF"].dev, 2)) + ")"  )

    report_dict[dataset].append( str(stat[1]["AU"].count) + " (" + str(round( 100.0 * float(stat[1]["AU"].count) / float(read_pairs), 2) ) + "%)" )
    report_dict[dataset].append( str(stat[1]["SP"].count) + " (" + str(round( 100.0 * float(stat[1]["SP"].count) / float(read_pairs), 2) ) + "%)" )

# total report
import report_maker
report_maker.do(report_dict, os.path.join(output_dir, 'report.horizontal'), os.path.join(output_dir, 'report'))

# clearing temp folder
shutil.rmtree(tmp_folder)

print("Done.")

