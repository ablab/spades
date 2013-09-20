#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import stat
import sys
import logging
import glob
import re
import options_storage

# constants to print and detect warnings and errors in logs
SPADES_PY_ERROR_MESSAGE = "== Error == "
SPADES_PY_WARN_MESSAGE = "== Warning == "
SPADES_ERROR_MESSAGE = " ERROR "
SPADES_WARN_MESSAGE = " WARN "


def error(err_str, log=None, prefix=SPADES_PY_ERROR_MESSAGE):
    if log:
        log.info("\n\n" + prefix + " " + err_str)
        log_warnings(log)
        log.info("\nIn case you have troubles running SPAdes, you can write to spades.support@bioinf.spbau.ru")
        log.info("Please provide us with params.txt and spades.log files from the output directory.")
    else:
        sys.stderr.write("\n\n" + prefix + " " + err_str + "\n\n")
        sys.stderr.write("In case you have troubles running SPAdes, you can write to spades.support@bioinf.spbau.ru\n")
        sys.stderr.write("Please provide us with params.txt and spades.log files from the output directory.\n")
        sys.stderr.flush()
    sys.exit(1)


def warning(warn_str, log=None, prefix="== Warning == "):
    if log:
        log.info("\n\n" + prefix + " " + warn_str + "\n\n")
    else:
        sys.stdout.write("\n\n" + prefix + " " + warn_str + "\n\n\n")
        sys.stdout.flush()


def check_python_version():
    if sys.version[0:3] not in options_storage.SUPPORTED_PYTHON_VERSIONS:
        error("python version " + sys.version[0:3] + " is not supported!\n" + \
              "Supported versions are " + ", ".join(options_storage.SUPPORTED_PYTHON_VERSIONS))


def check_file_existence(filename, message="", log=None):
    if not os.path.isfile(filename):
        error("file not found: %s (%s)" % (filename, message), log)
    return filename


def check_files_duplication(filenames, log):
    for filename in filenames:
        if filenames.count(filename) != 1:
            error("file %s was specified at least twice" % filename, log)


def check_reads_file_format(filename, message, only_assembler, log):
    if filename in options_storage.dict_of_prefixes:
        ext = options_storage.dict_of_prefixes[filename]
    else:
        ext = os.path.splitext(filename)[1]
        if ext.lower() == '.gz':
            ext = os.path.splitext(filename[:-len(ext)])[1] + ext
    if ext.lower() not in options_storage.ALLOWED_READS_EXTENSIONS:
        error("file with reads has unsupported format (only " + ", ".join(options_storage.ALLOWED_READS_EXTENSIONS) +
              " are supported): %s (%s)" % (filename, message), log)
    if not only_assembler and ext.lower() not in options_storage.BH_ALLOWED_READS_EXTENSIONS:
        error("to run read error correction, reads should be in FASTQ format (" +
              ", ".join(options_storage.BH_ALLOWED_READS_EXTENSIONS) +
              " are supported): %s (%s)" % (filename, message), log)


# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def process_subprocess_output(line):
    if sys.version.startswith('2.'):
        return line.rstrip()
    else: # sys.version.startswith('3.'):
        return str(line.rstrip(), 'utf-8')


def sys_call(cmd, log=None, cwd=None):
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    output = ''
    while not proc.poll():
        line = process_subprocess_output(proc.stdout.readline())
        if line:
            if log:
                log.info(line)
            else:
                output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = process_subprocess_output(line)
        if line:
            if log:
                log.info(line)
            else:
                output += line + "\n"

    if proc.returncode:
        error('system call for: "%s" finished abnormally, err code: %d' % (cmd, proc.returncode), log)
    return output


def universal_sys_call(cmd, log, out_filename=None, err_filename=None, cwd=None):
    '''
    Runs cmd and redirects stdout to out_filename (if specified), stderr to err_filename (if specified), or to log otherwise
    '''
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
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
                line = process_subprocess_output(proc.stdout.readline())
                if line:
                    log.info(line)
            if not err_filename:
                line = process_subprocess_output(proc.stderr.readline())
                if line:
                    log.info(line)
            if proc.returncode is not None:
                break

        if not out_filename:
            for line in proc.stdout.readlines():
                if line != '':
                    log.info(process_subprocess_output(line))
        if not err_filename:
            for line in proc.stderr.readlines():
                if line != '':
                    log.info(process_subprocess_output(line))
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


def get_warnings(log_filename):
    spades_py_warns = []
    spades_warns = []
    for line in open(log_filename, 'r'):
        if line.find(SPADES_PY_WARN_MESSAGE) != -1:
            line = line.replace(SPADES_PY_WARN_MESSAGE, '')
            spades_py_warns.append(' * ' + line.strip())
        elif line.find(SPADES_WARN_MESSAGE) != -1:
            spades_warns.append(' * ' + line.strip())
    return spades_py_warns, spades_warns


def get_logger_filename(log):
    log_file = None
    for h in log.__dict__['handlers']:
        if h.__class__.__name__ == 'FileHandler':
            log_file = h.baseFilename
    return log_file


def log_warnings(log):
    log_file = get_logger_filename(log)
    if not log_file:
        return False
    spades_py_warns, spades_warns = get_warnings(log_file)
    if spades_py_warns or spades_warns:
        log.info("\n======= SPAdes pipeline finished WITH WARNINGS!")
        warnings_filename = os.path.join(os.path.dirname(log_file), "warnings.log")
        warnings_handler = logging.FileHandler(warnings_filename, mode='w')
        log.addHandler(warnings_handler)
        #log.info("===== Warnings occurred during SPAdes run =====")
        log.info("")
        if spades_py_warns:
            log.info("=== Pipeline warnings:")
            for line in spades_py_warns:
                log.info(line)
        if spades_warns:
            log.info("=== Error correction and assembling warnings:")
            for line in spades_warns:
                log.info(line)
        log.info("======= Warnings saved to " + warnings_filename)
        log.removeHandler(warnings_handler)
        return True
    return False


def get_latest_dir(pattern):
    def atoi(text):
        if text.isdigit():
            return int(text)
        return text

    def natural_keys(text):
        return [atoi(c) for c in re.split('(\d+)', text)]

    latest_dir = None
    for dir_to_test in sorted(glob.glob(pattern), key=natural_keys, reverse=True):
        if os.path.isdir(dir_to_test):
            latest_dir = dir_to_test
            break
    return latest_dir


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
        if data.find(':') != -1 and ('.' + data[:data.find(':')]) in options_storage.ALLOWED_READS_EXTENSIONS:
            prefix = '.' + data[:data.find(':')]
            data = data[data.find(':') + 1:]
            options_storage.dict_of_prefixes[data] = prefix
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


def relative2abs_paths(dataset_data, dir_name):
    dir_name = os.path.abspath(dir_name)
    abs_paths_dataset_data = []
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith('reads'):
                abs_paths_reads = []
                for reads_file in value:
                    abs_path = os.path.join(dir_name, reads_file)
                    if reads_file in options_storage.dict_of_prefixes and abs_path != reads_file:
                        options_storage.dict_of_prefixes[abs_path] = options_storage.dict_of_prefixes[reads_file]
                        del options_storage.dict_of_prefixes[reads_file]
                    abs_paths_reads.append(abs_path)
                reads_library[key] = abs_paths_reads
        abs_paths_dataset_data.append(reads_library)
    return abs_paths_dataset_data


def check_dataset_reads(dataset_data, only_assembler, log):
    all_files = []
    for id, reads_library in enumerate(dataset_data):
        left_number = 0
        right_number = 0
        for key, value in reads_library.items():
            if key.endswith('reads'):
                for reads_file in value:
                    check_file_existence(reads_file, key + ', library number: ' + str(id + 1) +
                                         ', library type: ' + reads_library['type'], log)
                    check_reads_file_format(reads_file, key + ', library number: ' + str(id + 1) +
                                            ', library type: ' + reads_library['type'], only_assembler, log)
                    all_files.append(reads_file)
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


def split_interlaced_reads(dataset_data, dst, log):
    new_dataset_data = list()
    for reads_library in dataset_data:
        new_reads_library = dict(reads_library)
        for key, value in reads_library.items():
            if key == 'interlaced reads':
                if 'left reads' not in new_reads_library:
                    new_reads_library['left reads'] = []
                    new_reads_library['right reads'] = []
                for interlaced_reads in value:
                    if interlaced_reads in options_storage.dict_of_prefixes:
                        ext = options_storage.dict_of_prefixes[interlaced_reads]
                    else:
                        ext = os.path.splitext(interlaced_reads)[1]
                    was_compressed = False
                    if ext.endswith('.gz'):
                        was_compressed = True
                        import gzip
                        input_file = gzip.open(interlaced_reads, 'r')
                        ungzipped = os.path.splitext(interlaced_reads)[0]
                        out_basename, ext = os.path.splitext(os.path.basename(ungzipped))
                    else:
                        input_file = open(interlaced_reads, 'r')
                        out_basename, ext = os.path.splitext(os.path.basename(interlaced_reads))

                    if interlaced_reads in options_storage.dict_of_prefixes:
                        ext = options_storage.dict_of_prefixes[interlaced_reads]
                    if ext.lower().startswith('.fq') or ext.lower().startswith('.fastq'):
                        is_fasta_format = False
                        ext = '.fastq'
                    else:
                        is_fasta_format = True
                        ext = '.fasta'

                    out_left_filename = os.path.join(dst, out_basename + "_1" + ext)
                    out_right_filename = os.path.join(dst, out_basename + "_2" + ext)

                    if not (options_storage.continue_mode and os.path.isfile(out_left_filename) and os.path.isfile(out_right_filename)):
                        options_storage.continue_mode = False
                        log.info("== Splitting " + interlaced_reads + " into left and right reads (in " + dst + " directory)")
                        out_left_file = open(out_left_filename, 'w')
                        out_right_file = open(out_right_filename, 'w')
                        if sys.version.startswith('3.') and was_compressed:
                            for id, line in enumerate(input_file):
                                if (is_fasta_format and (id % 4 < 2)) or (not is_fasta_format and (id % 8 < 4)):
                                    out_left_file.write(str(line, 'utf-8'))
                                else:
                                    out_right_file.write(str(line, 'utf-8'))
                        else:
                            for id, line in enumerate(input_file):
                                if (is_fasta_format and (id % 4 < 2)) or (not is_fasta_format and (id % 8 < 4)):
                                    out_left_file.write(line)
                                else:
                                    out_right_file.write(line)
                        out_left_file.close()
                        out_right_file.close()

                    input_file.close()
                    new_reads_library['left reads'].append(out_left_filename)
                    new_reads_library['right reads'].append(out_right_filename)
                    if interlaced_reads in options_storage.dict_of_prefixes:
                        del options_storage.dict_of_prefixes[interlaced_reads]
                del new_reads_library['interlaced reads']
        new_dataset_data.append(new_reads_library)
    return new_dataset_data


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
                value = str(reads_library[reads_type])
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
        for i in range(0, len(seq), 60):
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