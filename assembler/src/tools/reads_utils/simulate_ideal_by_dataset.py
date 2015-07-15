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
sys.path.append(os.path.join(os.path.abspath(sys.path[0]), '..'))

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
skip_trimming = False

###################################################################

long_options = "output-dir= reference= thread-num= bin-size= kmer-size= skip-trimming".split()
short_options = "o:r:t:b:k:s"

def usage():
    print 'Estimation reads quality'
    print 'Usage:', sys.argv[0], ' [options described below] datasets description-file(s)'
    print ""
    print "Options with parameters:"
    print "-r\t--reference\tFile with reference genome (Mandatory parameter)"
    print "-o\t--output-dir\tDirectory to store all result files"
    print "-t\t--thread-num\tMax number of threads (default is " + str(thread_num) + ")"
    print "-s\t--skip-trimming\tSkip N-trimming for speed-up"
    
def check_file(f):
    if not os.path.isfile(f):
        print "Error - file not found:", f
        sys.exit(2)
    return f

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
    elif opt in ('-b', "--bin-size"):
        bin_size = int(arg)
        if bin_size < 1:
            bin_size = 1   
    elif opt in ('-k', "--kmer-size"):
        kmer = int(arg)
        if kmer < 1:
            kmer = 1     
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
print("Aligning")
report_dict["header"] += ["Total reads", "Aligned reads", "Unaligned reads"]
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = open(os.path.join(output_dir, dataset + ".log"),'w')
    align_err = open(os.path.join(output_dir, dataset + ".err"),'w') 
    reads_string = reduce(lambda x, y: x + ',' + y, datasets_dict[dataset])   
    subprocess.call([bowtie, '-c', '-q', '--suppress', '6,7,8', index, '-p', str(thread_num), reads_string], stdout=align_log, stderr=align_err)
    align_log.close()
    align_err.close() 

    align_err = open(os.path.join(output_dir, dataset + ".err"),'r') 
    for line in align_err:
        if line.startswith("# reads processed") or line.startswith("# reads with at least one") or line.startswith("# reads that failed"):
            report_dict[dataset].append( (line.split(':')[1]).strip() )
    align_err.close()       

# raw-single    
print("Parsing Bowtie log")
import raw_single
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = os.path.join(output_dir, dataset + ".log")
    raw_file  = os.path.join(output_dir, dataset + ".raw")
    raw_single.raw_single(align_log, raw_file)

# is form logs    
print("Retaining insert size")
import is_from_single_log
stat_dict = {}
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")
    align_log = os.path.join(output_dir, dataset + ".log")
    stat_dict[dataset] = is_from_single_log.dist_from_log(align_log, 1000)

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

print("Simulating ideal reads")
import fastaparser
import simulate_ideal_by_fasta
for dataset in datasets_dict.iterkeys():
    print("  " + dataset + "...")

    #print(str(gaps_dict[dataset][0][0][0]) + " " + str(gaps_dict[dataset][0][-1][1]) + " " + str(ref_len - 1) )

    circular = gaps_dict[dataset][0][0][0] <= 1 and gaps_dict[dataset][0][-1][1] >= ref_len - 1
    ins = int(stat_dict[dataset][1] + stat_dict[dataset][2] / 4.0)
    rl = stat_dict[dataset][0]
    
    #print (str(stat_dict[dataset][1]) + " " + str(stat_dict[dataset][2]))
    print ("     Read length: " + str(rl) + "; distance: " + str(ins) + "; circular: " + str(circular))

    chunks = gaps_dict[dataset][0]
    align_log = os.path.join(output_dir, dataset + ".log")
    last_gap = (chunks[-1][1], chunks[0][0])

    if circular or last_gap[1] + ref_len - last_gap[0] >= ins:
        last_gap = (0, 0)

    trusted_gaps = simulate_ideal_by_fasta.trusted_gaps(align_log, gaps_dict[dataset][1], last_gap, ins)

    if last_gap not in trusted_gaps:
        last_gap = (0,0)

    total_gaps = len(chunks) - 1
    if not circular:
        total_gaps += 1
    print ("     Total gaps between aligned regions: " + str(total_gaps) + "; gaps, that have read pairs spanning over: " + str(len(trusted_gaps)))

    trusted_chunks_index = [(0,0) for i in range(0, len(chunks))]
    for i in range(0, len(chunks) - 1):
        if (chunks[i][1], chunks[i + 1][0]) in trusted_gaps:
            trusted_chunks_index[i] = (chunks[i][1], chunks[i + 1][0])

    chunks_file  = os.path.join(output_dir, os.path.splitext(os.path.basename(reference))[0] + "gaps_" + dataset + ".fasta")
    fasta = fastaparser.read_fasta(chunks_file)

    singlef  = os.path.join(output_dir, dataset + "_single_reads.fasta")
    simulate_ideal_by_fasta.simulate_single(singlef, fasta, rl, circular)

    if ins < rl:
        continue

    pairedf  = os.path.join(output_dir, dataset + "_paired_reads.fasta")
    simulate_ideal_by_fasta.simulate_paired(pairedf, fasta, ins, rl, circular)

    gapf  = os.path.join(output_dir, dataset + "_gapped_reads.fasta")
    simulate_ideal_by_fasta.simulate_paired_over_gaps(gapf, fasta, chunks, trusted_chunks_index, last_gap, ref_len, ins, rl, circular)


# total report
import report_maker
report_maker.do(report_dict, os.path.join(output_dir, 'report.horizontal'), os.path.join(output_dir, 'report'))

# clearing temp folder
shutil.rmtree(tmp_folder)

print("Done.")

