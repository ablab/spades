#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import re
import support


def dataset_print(dataset):
    result = ""
    dataset_dict = dict(dataset)
    for key, value in dataset_dict.items():
        result += key + "\t" + value + "\n"
    return result


def generate_dataset(cfg):
    import process_cfg
    dataset_cfg = dict()
    dataset_cfg["single_cell"] = process_cfg.bool_to_str(cfg.single_cell)
    dataset_cfg["reads"] = os.path.join(cfg.output_dir, "corrected.yaml")
    if "reference_genome" in cfg.__dict__:
        dataset_cfg["reference_genome"] = cfg.reference_genome

    return dataset_print(dataset_cfg)


#### auxiliary function to manage input files 

def split_paired_file(input_filename, output_folder, log):
    ext = os.path.splitext(input_filename)[1]

    input_file = file
    out_basename = ""

    if ext == '.gz':
        import gzip

        input_file = gzip.open(input_filename, 'r')
        ungzipped = os.path.splitext(input_filename)[0]
        out_basename = os.path.splitext(os.path.basename(ungzipped))[0]
    else:
        input_file = open(input_filename, 'r')
        out_basename = os.path.splitext(os.path.basename(input_filename))[0]

    out_left_filename = os.path.join(output_folder, out_basename + "_1.fastq")
    out_right_filename = os.path.join(output_folder, out_basename + "_2.fastq")

    log.info("== Splitting " + input_filename + " into left and right reads")

    out_left_file = open(out_left_filename, 'w')
    out_right_file = open(out_right_filename, 'w')
    for id, line in enumerate(input_file):
        if id % 8 < 4:
            out_left_file.write(line)
        else:
            out_right_file.write(line)

    out_left_file.close()
    out_right_file.close()
    input_file.close()

    return [out_left_filename, out_right_filename]


def merge_paired_files(src_paired_reads, dst_paired_reads, output_folder, log):
    merged = []

    for i in [0, 1]:
        dst_basename = generate_paired_basename(dst_paired_reads[i], src_paired_reads[i], i)
        dst_filename = ""
        if dst_paired_reads[i].startswith(output_folder):
            dst_filename = os.path.join(os.path.dirname(dst_paired_reads[i]), dst_basename)
            os.rename(dst_paired_reads[i], dst_filename)
        else:
            import shutil

            dst_filename = os.path.join(output_folder, dst_basename)
            shutil.copy(dst_paired_reads[i], dst_filename)

        merged.append(dst_filename)

        log.info("== Merging " + src_paired_reads[i] + " and " + dst_paired_reads[
                                                              i] + " into one file with paired reads: " + dst_filename)

        src_file = open(src_paired_reads[i], 'r')
        dst_file = open(dst_filename, "a")
        dst_file.write(src_file.read())
        dst_file.close()
        src_file.close()

    return merged


def merge_single_files(src_single_read, dst_single_read, output_folder, log):
    dst_filename = ""
    if dst_single_read.startswith(output_folder):
        dst_filename = dst_single_read
    else:
        import shutil

        shutil.copy(dst_single_read, output_folder)
        dst_filename = os.path.join(output_folder, os.path.basename(dst_single_read))
    merged_filename = os.path.join(os.path.dirname(dst_filename),
        generate_unpaired_basename(src_single_read, dst_filename))

    log.info("== Merging " + src_single_read + " and " + dst_single_read + \
             " into one file with unpaired reads: " + merged_filename)

    src_file = open(src_single_read, 'r')
    dst_file = open(dst_filename, "a")
    dst_file.write(src_file.read())
    dst_file.close()
    src_file.close()
    os.rename(dst_filename, merged_filename)

    return merged_filename


def ungzip_if_needed(filename, output_folder, log):
    file_basename, file_extension = os.path.splitext(filename)
    if file_extension == ".gz":
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        ungzipped_filename = os.path.join(output_folder, os.path.basename(file_basename))
        ungzipped_file = open(ungzipped_filename, 'w')

        log.info("== Extracting " + filename + " into " + ungzipped_filename)

        import subprocess

        return_code = subprocess.call(['gunzip', filename, '-c'], stdout=ungzipped_file)
        support.verify(return_code == 0, log, "GZIP failed to extract " + filename + ". Maybe the archive is broken.")

        ungzipped_file.close()
        filename = ungzipped_filename

    return filename

####
