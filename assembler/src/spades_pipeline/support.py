#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import stat
import sys
import shutil


def verify(expr, log, message):
    if (not (expr)):
        log.info ("Assertion failed. Message: " + message)
        sys.exit(1)


def error(err_str, log=None, prefix="== Error == "):
    if log:
        log.info("\n\n" + prefix + " " + err_str + "\n")
        log.info("In case you have troubles running SPAdes, you can write to spades.support@bioinf.spbau.ru")
        log.info("Please provide us with params.txt and spades.log files from the output directory.\n")
    else:
        print >>sys.stderr, "\n\n" + prefix + " " + err_str + "\n"
        print >>sys.stderr, "In case you have troubles running SPAdes, you can write to spades.support@bioinf.spbau.ru"
        print >>sys.stderr, "Please provide us with params.txt and spades.log files from the output directory.\n"
    sys.exit(1)


def warning(warn_str, log=None, prefix="== Warning == "):
    if log:
        log.info("\n\n" + prefix + " " + warn_str + "\n\n")
    else:
        print "\n\n" + prefix + " " + warn_str + "\n\n"


def check_file_existence(filename, message="", log=None):
    if not os.path.isfile(filename):
        error("file not found: %s (%s)" % (filename, message), log)
    return filename


def check_files_duplication(filenames, log):
    for filename in filenames:
        if filenames.count(filename) != 1:
            error("file %s was specified at least twice" % filename, log)


def check_reads_file_format(filename, message, log):
    ext = os.path.splitext(filename)[1]
    if ext.lower() not in ['.fa', '.fasta', '.fq', '.fastq', '.gz']:
        error("file with reads has unsupported format (only .fa, .fasta, .fq,"
              " .fastq, .gz are supported): %s (%s)" % (filename, message), log)


def sys_call(cmd, log, cwd=None):
    import shlex
    import subprocess

    cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    while not proc.poll():
        line = proc.stdout.readline()
        if line != '':
            log.info(line.rstrip())
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        if line != '':
            log.info(line.rstrip())

    if proc.returncode:
        error('system call for: "%s" finished abnormally, err code: %d' % (cmd, proc.returncode), log)


def universal_sys_call(cmd, log, out_filename=None, err_filename=None, cwd=None):
    '''
    Runs cmd and redirects stdout to out_filename (if specified), stderr to err_filename (if specified), or to log otherwise
    '''
    import shlex
    import subprocess

    cmd_list = shlex.split(cmd)

    if out_filename:
        stdout = open(out_filename, 'w')
    else:
        stdout = subprocess.PIPE
    if err_filename:
        stderr = open(err_filename, 'w')
    else:
        stderr = subprocess.PIPE

    proc = subprocess.Popen(cmd_list, stdout=stdout, stderr=stderr, cwd=cwd)

    if log and (not out_filename or not err_filename):
        while not proc.poll():
            if not out_filename:
                line = proc.stdout.readline()
                if line != '':
                    log.info(line.rstrip())
            if not err_filename:
                line = proc.stderr.readline()
                if line != '':
                    log.info(line.rstrip())
            if proc.returncode is not None:
                break

        if not out_filename:
            for line in proc.stdout.readlines():
                if line != '':
                    log.info(line.rstrip())
        if not err_filename:
            for line in proc.stderr.readlines():
                if line != '':
                    log.info(line.rstrip())
    else:
        proc.wait()

    if out_filename:
        stdout.close()
    if err_filename:
        stderr.close()


def save_data_to_file(data, file):
    output = open(file, 'wb')
    output.write(data.read())
    output.close()
    os.chmod(file, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


### START for processing YAML files
def get_lib_type_and_number(option):
    lib_type = 'pe'
    lib_number = 1
    if option.startswith('--mp'):
        lib_type = 'mp'
    if option.startswith('--mp') or option.startswith('--pe'): # don't process simple -1, -2, -s, --12 options
        lib_number = int(option[4])
    return lib_type, lib_number


def get_data_type(option):
    if option.endswith('-12'):
        data_type = 'interlaced reads'
    elif option.endswith('-1'):
        data_type = 'left reads'
    elif option.endswith('-2'):
        data_type = 'right reads'
    elif option.endswith('-s'):
        data_type = 'single reads'
    else: # -rf, -ff, -fr
        data_type = 'orientation'
    return data_type


def add_to_dataset(option, data, dataset_data):
    lib_type, lib_number = get_lib_type_and_number(option)
    data_type = get_data_type(option)
    if data_type == 'orientation':
        data = option[-2:]
    total_libs_number = len(dataset_data) / 2

    if lib_type == 'pe':
        record_id = lib_number - 1
    else: # mate-pairs
        record_id = total_libs_number + lib_number - 1

    if not dataset_data[record_id]: # setting default values for a new record
        if lib_type == 'pe':
            dataset_data[record_id]['type'] = 'paired-end'
        else:
            dataset_data[record_id]['type'] = 'mate-pairs'
    if data_type.endswith('reads'): # reads are stored as lists
        if data_type in dataset_data[record_id]:
            dataset_data[record_id][data_type].append(data)
        else:
            dataset_data[record_id][data_type] = [data]
    else: # other values are stored as plain strings
        dataset_data[record_id][data_type] = data


def correct_dataset(dataset_data):
    # removing empty reads libraries
    corrected_dataset_data = []
    for reads_library in dataset_data:
        if not reads_library:
            continue
        has_reads = False
        has_paired_reads = False
        for key in reads_library.keys():
            if key.endswith('reads'):
                has_reads = True
            if key in ['interlaced reads', 'left reads', 'right reads']:
                has_paired_reads = True
                break
        if not has_reads:
            continue
        if not has_paired_reads and reads_library['type'] == 'paired-end':
            reads_library['type'] = 'single'
            if 'orientation' in reads_library:
                del reads_library['orientation']
        if 'orientation' not in reads_library:
            if reads_library['type'] == 'paired-end':
                reads_library['orientation'] = 'fr'
            elif reads_library['type'] == 'mate-pairs':
                reads_library['orientation'] = 'rf'
        corrected_dataset_data.append(reads_library)
    return corrected_dataset_data


def check_dataset_reads(dataset_data, log):
    all_files = []
    for id, reads_library in enumerate(dataset_data):
        left_number = 0
        right_number = 0
        for key, value in reads_library.items():
            if key.endswith('reads'):
                for reads_file in value:
                    check_file_existence(os.path.abspath(reads_file), key + ', library number: ' + str(id + 1) +
                                         ', library type: ' + reads_library['type'], log)
                    check_reads_file_format(os.path.abspath(reads_file), key + ', library number: ' + str(id + 1) +
                                            ', library type: ' + reads_library['type'], log)
                    all_files.append(os.path.abspath(reads_file))
                if key == 'left reads':
                    left_number = len(value)
                elif key == 'right reads':
                    right_number = len(value)
        if left_number != right_number:
            error('the number of files with left paired reads is not equal to the number of files '
                  'with right paired reads (library number: ' + str(id + 1) +
                  ', library type: ' + reads_library['type'] + ')!', log)
    if not len(all_files):
        error("You should specify at least one file with reads!", log)
    check_files_duplication(all_files, log)


def dataset_has_only_mate_pairs_libraries(dataset_data):
    for reads_library in dataset_data:
        if reads_library['type'] != 'mate-pairs':
            return False
    return True


def dataset_has_paired_reads(dataset_data):
    for reads_library in dataset_data:
        if reads_library['type'] in ['paired-end', 'mate-pairs']:
            return True
    return False


def dataset_has_interlaced_reads(dataset_data):
    for reads_library in dataset_data:
        if 'interlaced reads' in reads_library:
            return True
    return False


def move_dataset_files(dataset_data, dst, log, gzip=False):
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith('reads'):
                moved_reads_files = []
                for reads_file in value:
                    dst_filename = os.path.join(dst, os.path.basename(reads_file))
                    # TODO: fix problem with files with the same basenames in Hammer binary!
                    if not os.path.isfile(reads_file):
                        if (not gzip and os.path.isfile(dst_filename)) or (gzip and os.path.isfile(dst_filename + '.gz')):
                            warning('file with corrected reads (' + reads_file + ') is the same in several libraries')
                            if gzip:
                                dst_filename += '.gz'
                        else:
                            error('something went wrong and file with corrected reads (' + reads_file + ') is missing!')
                    else:
                        shutil.move(reads_file, dst_filename)
                        if gzip:
                            #log.info('Compressing ' + dst_filename + ' into ' + dst_filename + '.gz')
                            sys_call('gzip -f -9 ' + dst_filename, log)
                            dst_filename += '.gz'
                    moved_reads_files.append(dst_filename)
                reads_library[key] = moved_reads_files


def split_interlaced_reads(dataset_data, dst, log):
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key == 'interlaced reads':
                if 'left reads' not in reads_library:
                    reads_library['left reads'] = []
                    reads_library['right reads'] = []
                for interlaced_reads in value:
                    ext = os.path.splitext(interlaced_reads)[1]
                    if ext == '.gz':
                        import gzip
                        input_file = gzip.open(interlaced_reads, 'r')
                        ungzipped = os.path.splitext(interlaced_reads)[0]
                        out_basename = os.path.splitext(os.path.basename(ungzipped))[0]
                    else:
                        input_file = open(interlaced_reads, 'r')
                        out_basename = os.path.splitext(os.path.basename(interlaced_reads))[0]
                    out_left_filename = os.path.join(dst, out_basename + "_1.fastq")
                    out_right_filename = os.path.join(dst, out_basename + "_2.fastq")

                    log.info("== Splitting " + interlaced_reads + " into left and right reads (in " + dst + " directory)")
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
                    reads_library['left reads'].append(out_left_filename)
                    reads_library['right reads'].append(out_right_filename)
                del reads_library['interlaced reads']


def pretty_print_reads(dataset_data, log, indent='    '):
    READS_TYPES = ['left reads', 'right reads', 'interlaced reads', 'single reads']
    for id, reads_library in enumerate(dataset_data):
        log.info(indent + 'Library number: ' + str(id + 1) + ', library type: ' + reads_library['type'])
        if 'orientation' in reads_library:
            log.info(indent + '  orientation: ' + reads_library['orientation'])
        for reads_type in READS_TYPES:
            if reads_type not in reads_library:
                value = 'not specified'
            else:
                value = str(map(os.path.abspath, reads_library[reads_type]))
            log.info(indent + '  ' + reads_type + ': ' + value)
### END: for processing YAML files


def read_fasta(filename):
    """
        Returns list of FASTA entries (in tuples: name, seq)
    """
    res_name = []
    res_seq = []
    first = True
    seq = ''

    for line in open(filename):
        if line[0] == '>':
            res_name.append(line.strip())
            if not first:
                res_seq.append(seq)
            else:
                first = False
            seq = ''
        else:
            seq += line.strip()
    res_seq.append(seq)
    return zip(res_name, res_seq)


def write_fasta(filename, fasta):
    outfile = open(filename, 'w')
    for name, seq in fasta:
        outfile.write(name + '\n')
        for i in xrange(0,len(seq),60):
            outfile.write(seq[i : i + 60] + '\n')
    outfile.close()


def break_scaffolds(input_filename, threshold, output_filename):
    new_fasta = []
    for id, (name, seq) in enumerate(read_fasta(input_filename)):
        i = 0
        cur_contig_number = 1
        cur_contig_start = 0
        while (i < len(seq)) and (seq.find("N", i) != -1):
            start = seq.find("N", i)
            end = start + 1
            while (end != len(seq)) and (seq[end] == 'N'):
                end += 1

            i = end + 1
            if (end - start) >= threshold:
                new_fasta.append((name.split()[0] + "_" + str(cur_contig_number), seq[cur_contig_start:start]))
                cur_contig_number += 1
                cur_contig_start = end

        new_fasta.append((name.split()[0] + "_" + str(cur_contig_number), seq[cur_contig_start:]))

    write_fasta(output_filename, new_fasta)