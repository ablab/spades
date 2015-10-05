############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import support
import os
import shutil

def align_bwa_pe_lib(command, index, reads_file1, reads_file2, work_dir, log, threads = 1):
    log.info("Aligning paired-end library")
    log.info("Left reads: " + reads_file1)
    log.info("Right reads: " + reads_file2)
    log.info("Output directory: " + work_dir)
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    log_file = os.path.join(work_dir, "output.log")
    err_log_file = os.path.join(work_dir, "err_output")
    result = os.path.join(work_dir, "alignment.sam")
    log.info("Starting alignment of reads using bwa. See detailed log in " + log_file)
    log.info("Starting read alignment. See detailed log in " + log_file)
    support.universal_sys_call([command, "mem", "-t", str(threads), "-S", "-M", index, reads_file1, reads_file2], log, result, err_log_file)
    log.info("Done. See result in " + result)
    return result


def index_bwa(command, log, reference, work_dir, algorithm = "is"):
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    log.info("Constructing bwa index")
    index = os.path.join(work_dir, "index")
    log_file = os.path.join(work_dir, "output.log")
    err_log_file = os.path.join(work_dir, "err_output")
    log.info(" ".join([command, "index", "-a", algorithm, "-p", index, reference]))
    support.universal_sys_call([command, "index", "-a", algorithm, "-p", index, reference], log, log_file, err_log_file)
    log.info("Index constructed.")
    return index


def align_bwa_pe_libs(command, index, reads, work_dir, log, threads):
    log.info("===== Starting read alignment")
    result = []
    lib_num = 1
    for left_reads, right_reads in reads:
        lib_dir = os.path.join(work_dir, str(lib_num))
        result.append(
            align_bwa_pe_lib(command, index, left_reads, right_reads, lib_dir, log, threads))
        lib_num += 1
    log.info("===== Read alignment finished. See result in " + work_dir)
    return result



def align_bwa(command, reference, dataset, work_dir, log = None, threads = 1):
    if log == None:
        log = logging.getLogger('')
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    index_dir = os.path.join(work_dir, "index")
    index = index_bwa(command, log, reference, index_dir)
    reads = [(lib["left reads"][0], lib["right reads"][0]) for lib in dataset]
    return align_bwa_pe_libs(command, index, reads, work_dir, log, threads)
