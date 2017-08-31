#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2016 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import stat
import sys
import logging
import glob
import re
import gzip
import tempfile
import shutil
import options_storage
import itertools
from common import SeqIO
import math
from os.path import abspath, expanduser, join
from distutils.version import LooseVersion

# constants to print and detect warnings and errors in logs
SPADES_PY_ERROR_MESSAGE = "== Error == "
SPADES_PY_WARN_MESSAGE = "== Warning == "
SPADES_ERROR_MESSAGE = " ERROR "
SPADES_WARN_MESSAGE = " WARN "
# for correct warnings detection in case of continue_mode
continue_logfile_offset = None
# for removing tmp_dir even if error occurs
current_tmp_dir = None


def error(err_str, log=None, dipspades=False, prefix=SPADES_PY_ERROR_MESSAGE):
    if not dipspades:
        binary_name = "SPAdes"
    else:
        binary_name = "dipSPAdes"
    if log:
        log.info("\n\n" + prefix + " " + err_str)
        log_warnings(log, with_error=True)
        log.info("\nIn case you have troubles running " + binary_name + ", you can write to spades.support@cab.spbu.ru")
        log.info("or report an issue on our GitHub repository github.com/ablab/spades")
        log.info("Please provide us with params.txt and " + binary_name.lower() + ".log files from the output directory.")
    else:
        sys.stderr.write("\n\n" + prefix + " " + err_str + "\n\n")
        sys.stderr.write("\nIn case you have troubles running " + binary_name + ", you can write to spades.support@cab.spbu.ru\n")
        sys.stderr.write("or report an issue on our GitHub repository github.com/ablab/spades\n")
        sys.stderr.write("Please provide us with params.txt and " + binary_name.lower() + ".log files from the output directory.\n")
        sys.stderr.flush()
    if current_tmp_dir and os.path.isdir(current_tmp_dir):
        shutil.rmtree(current_tmp_dir)
    sys.exit(1)


def warning(warn_str, log=None, prefix="== Warning == "):
    if log:
        log.info("\n\n" + prefix + " " + warn_str + "\n\n")
    else:
        sys.stdout.write("\n\n" + prefix + " " + warn_str + "\n\n\n")
        sys.stdout.flush()


def check_python_version():
    def __next_version(version):
        components = version.split('.')
        for i in reversed(range(len(components))):
            if components[i].isdigit():
                components[i] = str(int(components[i]) + 1)
                break
        return '.'.join(components)

    current_version = sys.version.split()[0]
    supported_versions_msg = []
    for supported_versions in options_storage.SUPPORTED_PYTHON_VERSIONS:
        major = supported_versions[0]
        if '-' in supported_versions:  # range
            min_inc, max_inc = supported_versions.split('-')
        elif supported_versions.endswith('+'):  # half open range
            min_inc, max_inc = supported_versions[:-1], major
        else:  # exact version
            min_inc = max_inc = supported_versions
        max_exc = __next_version(max_inc)
        supported_versions_msg.append("Python%s: %s" % (major, supported_versions.replace('+', " and higher")))
        if LooseVersion(min_inc) <= LooseVersion(current_version) < LooseVersion(max_exc):
            return True
    error("Python version " + current_version + " is not supported!\n" +
          "Supported versions are " + ", ".join(supported_versions_msg))


def get_spades_binaries_info_message():
    return "You can obtain SPAdes binaries in one of two ways:" +\
           "\n1. Download them from http://cab.spbu.ru/software/spades/" +\
           "\n2. Build source code with ./spades_compile.sh script"


def check_binaries(binary_dir, log):
    for binary in ["hammer", "ionhammer", "spades", "bwa-spades", "dipspades"]:
        binary_path = os.path.join(binary_dir, binary)
        if not os.path.isfile(binary_path):
            error("SPAdes binaries not found: " + binary_path + "\n" + get_spades_binaries_info_message(), log)


def check_file_existence(input_filename, message="", log=None, dipspades=False):
    filename = abspath(expanduser(input_filename))
    check_path_is_ascii(filename, message)
    if not os.path.isfile(filename):
        error("file not found: %s (%s)" % (filename, message), log=log, dipspades=dipspades)
    options_storage.dict_of_rel2abs[input_filename] = filename
    return filename


def check_dir_existence(input_dirname, message="", log=None, dipspades=False):
    dirname = abspath(expanduser(input_dirname))
    check_path_is_ascii(dirname, message)
    if not os.path.isdir(dirname):
        error("directory not found: %s (%s)" % (dirname, message), log=log, dipspades=dipspades)
    options_storage.dict_of_rel2abs[input_dirname] = dirname
    return dirname


def check_path_is_ascii(path, message=""):
    if not is_ascii_string(path):
        error("path contains non-ASCII characters: %s (%s)" % (path, message))


def ensure_dir_existence(dirname):
    if os.path.isfile(dirname):
        os.remove(dirname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def recreate_dir(dirname):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)


def check_files_duplication(filenames, log):
    for filename in filenames:
        if filenames.count(filename) != 1:
            error("file %s was specified at least twice" % filename, log)


def check_reads_file_format(filename, message, only_assembler, library_type, log):
    if filename in options_storage.dict_of_prefixes:
        ext = options_storage.dict_of_prefixes[filename]
    else:
        ext = os.path.splitext(filename)[1]
        if ext.lower() == '.gz':
            pre_ext = os.path.splitext(filename[:-len(ext)])[1]
            if (pre_ext + ext).lower() in options_storage.ALLOWED_READS_EXTENSIONS:
                ext = pre_ext + ext
            else: # allows ".fastq.1.gz" like extensions
                pre_pre_ext = os.path.splitext(filename[:-len(pre_ext + ext)])[1]
                ext = pre_pre_ext + ext
    if ext.lower() not in options_storage.ALLOWED_READS_EXTENSIONS:
        error("file with reads has unsupported format (only " + ", ".join(options_storage.ALLOWED_READS_EXTENSIONS) +
              " are supported): %s (%s)" % (filename, message), log)
    if not only_assembler and ext.lower() not in options_storage.BH_ALLOWED_READS_EXTENSIONS and \
       library_type not in options_storage.LONG_READS_TYPES:
        error("to run read error correction, reads should be in FASTQ format (" +
              ", ".join(options_storage.BH_ALLOWED_READS_EXTENSIONS) +
              " are supported): %s (%s)" % (filename, message), log)
    if library_type.endswith("contigs") and ext.lower() not in options_storage.CONTIGS_ALLOWED_READS_EXTENSIONS:
        error("file with " + library_type + " should be in FASTA format  (" +
              ", ".join(options_storage.CONTIGS_ALLOWED_READS_EXTENSIONS) +
              " are supported): %s (%s)" % (filename, message), log)


# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    elif "PATH" in os.environ:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def get_available_memory():
    mem_info_filename = "/proc/meminfo"
    avail_mem_header = "MemTotal:"
    if os.path.isfile(mem_info_filename):
        try:
            for line in open(mem_info_filename):
                if line.startswith(avail_mem_header):
                    avail_mem = int(line[len(avail_mem_header):].split()[0]) # in kB
                    avail_mem /= 1024 * 1024 # in GB
                    return avail_mem
        except ValueError:
            return None
        except IOError:
            return None
    return None


# based on http://stackoverflow.com/questions/196345/how-to-check-if-a-string-in-python-is-in-ascii
def is_ascii_string(line):
    try:
        line.encode('ascii')
    except UnicodeDecodeError:  # python2
        return False
    except UnicodeEncodeError:  # python3
        return False
    else:
        return True


def process_readline(line, is_python3=sys.version.startswith('3.')):
    if is_python3:
        return str(line, 'utf-8').rstrip()
    return line.rstrip()


def process_spaces(str):
    if str.find(" ") != -1:
        str = '"' + str + '"'
    return str


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
        line = process_readline(proc.stdout.readline())
        if line:
            if log:
                log.info(line)
            else:
                output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = process_readline(line)
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
                line = process_readline(proc.stdout.readline())
                if line:
                    log.info(line)
            if not err_filename:
                line = process_readline(proc.stderr.readline())
                if line:
                    log.info(line)
            if proc.returncode is not None:
                break

        if not out_filename:
            for line in proc.stdout.readlines():
                if line != '':
                    log.info(process_readline(line))
        if not err_filename:
            for line in proc.stderr.readlines():
                if line != '':
                    log.info(process_readline(line))
    else:
        proc.wait()

    if out_filename:
        stdout.close()
    if err_filename:
        stderr.close()
    if proc.returncode:
        error('system call for: "%s" finished abnormally, err code: %d' % (cmd, proc.returncode), log)


def save_data_to_file(data, file):
    output = open(file, 'wb')
    output.write(data.read())
    output.close()
    os.chmod(file, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def get_important_messages_from_log(log_filename, warnings=True):
    def already_saved(list_to_check, suffix): # for excluding duplicates (--continue-from may cause them)
        for item in list_to_check:
            if item.endswith(suffix):
                return True
        return False

    if warnings:
        spades_py_message = SPADES_PY_WARN_MESSAGE
        spades_message = SPADES_WARN_MESSAGE
    else:  # errors
        spades_py_message = SPADES_PY_ERROR_MESSAGE
        spades_message = SPADES_ERROR_MESSAGE

    ### for capturing correct warnings in case of continue_mode
    if continue_logfile_offset:
        continued_log = open(log_filename, 'r')
        continued_log.seek(continue_logfile_offset)
        continued_stage_phrase = continued_log.readline()
        while not continued_stage_phrase.strip():
            continued_stage_phrase = continued_log.readline()
        lines_to_check = continued_log.readlines()
        continued_log.close()

        all_lines = open(log_filename, 'r').readlines()
        failed_stage_index = all_lines.index(continued_stage_phrase)
        lines_to_check = all_lines[:failed_stage_index] + lines_to_check
    else:
        lines_to_check = open(log_filename, 'r').readlines()

    spades_py_msgs = []
    spades_msgs = []
    IMPORTANT_MESSAGE_SUMMARY_PREFIX = ' * '
    for line in lines_to_check:
        if line.startswith(IMPORTANT_MESSAGE_SUMMARY_PREFIX):
            continue
        if line.find(spades_py_message) != -1:
            suffix = line[line.find(spades_py_message) + len(spades_py_message):].strip()
            line = line.replace(spades_py_message, '').strip()
            if not already_saved(spades_py_msgs, suffix):
                spades_py_msgs.append(IMPORTANT_MESSAGE_SUMMARY_PREFIX + line)
        elif line.find(spades_message) != -1:
            suffix = line[line.find(spades_message) + len(spades_message):].strip()
            line = line.strip()
            if not already_saved(spades_msgs, suffix):
                spades_msgs.append(IMPORTANT_MESSAGE_SUMMARY_PREFIX + line)
    return spades_py_msgs, spades_msgs


def get_logger_filename(log):
    log_file = None
    for h in log.__dict__['handlers']:
        if h.__class__.__name__ == 'FileHandler':
            log_file = h.baseFilename
    return log_file


def log_warnings(log, with_error=False):
    log_file = get_logger_filename(log)
    if not log_file:
        return False
    for h in log.__dict__['handlers']:
        h.flush()
    spades_py_warns, spades_warns = get_important_messages_from_log(log_file, warnings=True)
    if spades_py_warns or spades_warns:
        if with_error:
            log.info("\n======= SPAdes pipeline finished abnormally and WITH WARNINGS!")
        else:
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
        if with_error:
            spades_py_errors, spades_errors = get_important_messages_from_log(log_file, warnings=False)
            log.info("")
            log.info("=== ERRORs:")
            for line in (spades_errors + spades_py_errors):
                log.info(line)
        return True
    return False


def continue_from_here(log):
    if options_storage.continue_mode:
        options_storage.continue_mode = False
        log_filename = get_logger_filename(log)
        if log_filename:
            log_file = open(log_filename, 'r')
            log_file.seek(0, 2) # seek to the end of file
            global continue_logfile_offset
            continue_logfile_offset = log_file.tell()


def finish_here(log):
    log.info("\n======= Skipping the rest of SPAdes pipeline (--stop-after was set to '%s'). "
             "You can continue later with --continue or --restart-from options\n" % options_storage.stop_after)
    options_storage.run_completed = True


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


def get_tmp_dir(prefix="", base_dir=None):
    global current_tmp_dir

    if not base_dir:
        base_dir = options_storage.tmp_dir
    if not os.path.isdir(base_dir):
        os.makedirs(base_dir)
    current_tmp_dir = tempfile.mkdtemp(dir=base_dir, prefix=prefix)
    return current_tmp_dir


### START for processing YAML files
def get_short_reads_type(option):
    for short_reads_type in options_storage.SHORT_READS_TYPES.keys():
        if option.startswith('--' + short_reads_type):
            # additional check to except collisions with LONG_READS_TYPES, e.g. --s<#> and --sanger
            if option[len('--' + short_reads_type):len('--' + short_reads_type) + 1].isdigit():
                return short_reads_type
    return None


def get_long_reads_type(option):
    for long_reads_type in options_storage.LONG_READS_TYPES:
        if option.startswith('--') and option in ("--" + long_reads_type):
            return long_reads_type
    return None


def is_single_read_type(option):
    return option.startswith('--s') and option[3:].isdigit()


def get_lib_type_and_number(option):
    # defaults for simple -1, -2, -s, --12 options
    lib_type = 'pe'
    lib_number = 1

    if get_short_reads_type(option):
        lib_type = get_short_reads_type(option)
        lib_number = int(option[re.search("\d", option).start()])
    elif get_long_reads_type(option):
        lib_type = get_long_reads_type(option)
    return lib_type, lib_number


def get_data_type(option):
    if option.endswith('-12'):
        data_type = 'interlaced reads'
    elif option.endswith('-1'):
        data_type = 'left reads'
    elif option.endswith('-2'):
        data_type = 'right reads'
    elif option.endswith('-s') or is_single_read_type(option) or get_long_reads_type(option):
        data_type = 'single reads'
    else: # -rf, -ff, -fr
        data_type = 'orientation'
    return data_type


def get_option_prefix(data):
    prefix = None
    if data.find(':') != -1 and ('.' + data[:data.find(':')]) in options_storage.ALLOWED_READS_EXTENSIONS:
        prefix = data[:data.find(':')]
        data = data[data.find(':') + 1:]
    return data, prefix


def add_to_dataset(option, data, dataset_data):
    lib_type, lib_number = get_lib_type_and_number(option)
    data_type = get_data_type(option)
    if data_type == 'orientation':
        data = option[-2:]

    if lib_type in options_storage.SHORT_READS_TYPES:
        record_id = options_storage.MAX_LIBS_NUMBER * sorted(options_storage.SHORT_READS_TYPES.keys()).index(lib_type) \
                    + lib_number - 1
    elif lib_type in options_storage.LONG_READS_TYPES:
        record_id = options_storage.MAX_LIBS_NUMBER * len(options_storage.SHORT_READS_TYPES.keys()) \
                    + options_storage.LONG_READS_TYPES.index(lib_type)
    else:
        error("can't detect library type from option %s!" % option)

    if not dataset_data[record_id]:  # setting default values for a new record
        if lib_type in options_storage.SHORT_READS_TYPES:
            dataset_data[record_id]['type'] = options_storage.SHORT_READS_TYPES[lib_type]
        else:
            dataset_data[record_id]['type'] = lib_type
    if data_type.endswith('reads'):
        data, prefix = get_option_prefix(data)
        if prefix:
            options_storage.dict_of_prefixes[data] = '.' + prefix
        if data_type in dataset_data[record_id]:
            dataset_data[record_id][data_type].append(data)
        else:
            dataset_data[record_id][data_type] = [data]
    else:  # other values are stored as plain strings
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
            if reads_library['type'] == 'paired-end' or reads_library['type'] == 'hq-mate-pairs':
                reads_library['orientation'] = 'fr'
            elif reads_library['type'] == 'mate-pairs':
                reads_library['orientation'] = 'rf'
        corrected_dataset_data.append(reads_library)
    return corrected_dataset_data


def relative2abs_paths(dataset_data, dirname):
    dirname = abspath(expanduser(dirname))
    abs_paths_dataset_data = []
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith('reads'):
                abs_paths_reads = []
                for reads_file in value:
                    abs_path = abspath(join(dirname, expanduser(reads_file)))
                    options_storage.dict_of_rel2abs[reads_file] = abs_path
                    if reads_file in options_storage.dict_of_prefixes and abs_path != reads_file:
                        options_storage.dict_of_prefixes[abs_path] = options_storage.dict_of_prefixes[reads_file]
                        del options_storage.dict_of_prefixes[reads_file]
                    abs_paths_reads.append(abs_path)
                reads_library[key] = abs_paths_reads
        abs_paths_dataset_data.append(reads_library)
    return abs_paths_dataset_data


def get_reads_length(dataset_data, log, num_checked=10 ** 4, diff_len_allowable=25):
    max_reads_lenghts = [get_max_reads_length(reads_file, log, num_checked) for reads_file in get_reads_files(dataset_data)]

    avg_len = sum(max_reads_lenghts) / len(max_reads_lenghts)
    for max_len in max_reads_lenghts:
        if math.fabs(max_len - avg_len) > diff_len_allowable:
            warning('Read lengths differ more than allowable. Length: ' + str(max_len) + '. Avg. length: ' + str(avg_len) + '.', log)
    reads_length = min(max_reads_lenghts)
    log.info('Reads length: ' + str(reads_length))
    return reads_length


def get_reads_files(dataset_data):
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith('reads'):
                for reads_file in value:
                    yield reads_file


def get_max_reads_length(reads_file, log, num_checked):
    file_type = SeqIO.get_read_file_type(reads_file)
    if not file_type:
        error('Incorrect type of reads file: ' + reads_file, log)

    max_reads_length = max([len(rec) for rec in itertools.islice(SeqIO.parse(SeqIO.Open(reads_file, "r"), file_type), num_checked)])
    log.info('Max reads length: ' + str(max_reads_length))
    return max_reads_length


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
                                            ', library type: ' + reads_library['type'], only_assembler, reads_library['type'], log)
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


def check_single_reads_in_options(options, log):
    only_old_style_options = True
    old_style_single_reads = False
    for option in options:
        if option not in options_storage.reads_options:
            continue
        if option in options_storage.OLD_STYLE_READS_OPTIONS:
            if option == '-s':
                old_style_single_reads = True
        else:
            only_old_style_options = False
    if not only_old_style_options and old_style_single_reads:
        warning("It is recommended to specify single reads with --pe<#>-s, --mp<#>-s, --hqmp<#>-s, "
                "or --s<#> option instead of -s!", log)


def get_lib_ids_by_type(dataset_data, types):
    if type(types) is not list:
        types = [types]
    lib_ids = []
    for id, reads_library in enumerate(dataset_data):
        if reads_library['type'] in types:
            lib_ids.append(id)
    return lib_ids


def get_libs_by_type(dataset_data, types):
    ids = get_lib_ids_by_type(dataset_data, types)
    result = []
    for id in ids:
        result.append(dataset_data[id])
    return result


def rm_libs_by_type(dataset_data, types):
    ids = get_lib_ids_by_type(dataset_data, types)
    for id in sorted(ids, reverse=True):
        del dataset_data[id]
    return dataset_data


def dataset_is_empty(dataset_data):
    for reads_library in dataset_data:
        if reads_library:
            return False
    return True


def dataset_has_gzipped_reads(dataset_data):
    for reads_library in dataset_data:
        for key in reads_library:
            if key.endswith('reads'):
                for reads_file in reads_library[key]:
                    if reads_file.endswith('.gz'):
                        return True
    return False


def dataset_has_interlaced_reads(dataset_data):
    for reads_library in dataset_data:
        if 'interlaced reads' in reads_library:
            return True
    return False


def dataset_has_additional_contigs(dataset_data):
    for reads_library in dataset_data:
        if reads_library['type'].endswith('contigs'):
            return True
    return False


def dataset_has_nxmate_reads(dataset_data):
    for reads_library in dataset_data:
        if reads_library['type'] == 'nxmate':
            return True
    return False


def process_Ns_in_additional_contigs(dataset_data, dst, log):
    new_dataset_data = list()
    for reads_library in dataset_data:
        new_reads_library = dict(reads_library)
        if reads_library["type"].endswith("contigs"):
            new_entry = []
            for contigs in reads_library["single reads"]:
                if contigs in options_storage.dict_of_prefixes:
                    ext = options_storage.dict_of_prefixes[contigs]
                    basename = contigs
                else:
                    basename, ext = os.path.splitext(contigs)
                gzipped = False
                if ext.endswith('.gz'):
                    gzipped = True
                    if contigs not in options_storage.dict_of_prefixes:
                        basename, _ = os.path.splitext(basename)
                modified, new_fasta = break_scaffolds(contigs, options_storage.THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS,
                    replace_char='A', gzipped=gzipped)
                if modified:
                    if not os.path.isdir(dst):
                        os.makedirs(dst)
                    new_filename = os.path.join(dst, os.path.basename(basename) + '.fasta')
                    if contigs in options_storage.dict_of_prefixes:
                        del options_storage.dict_of_prefixes[contigs]
                    log.info("== Processing additional contigs (%s): changing Ns to As and "
                             "splitting by continues (>= %d) Ns fragments (results are in %s directory)" % (contigs,
                             options_storage.THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS, dst))
                    write_fasta(new_filename, new_fasta)
                    new_entry.append(new_filename)
                else:
                    new_entry.append(contigs)
            new_reads_library["single reads"] = new_entry
        new_dataset_data.append(new_reads_library)
    return new_dataset_data


def split_interlaced_reads(dataset_data, dst, log):
    def write_single_read(in_file, out_file, read_name=None, is_fastq=False, is_python3=False):
        if read_name is None:
            read_name = process_readline(in_file.readline(), is_python3)
        if not read_name:
            return ''  # no next read
        read_value = process_readline(in_file.readline(), is_python3)
        line = process_readline(in_file.readline(), is_python3)
        fpos = in_file.tell()
        while (is_fastq and not line.startswith('+')) or (not is_fastq and not line.startswith('>')):
            read_value += line
            line = process_readline(in_file.readline(), is_python3)
            if not line:
                if fpos == in_file.tell():
                    break
                fpos = in_file.tell()
        out_file.write(read_name + '\n')
        out_file.write(read_value + '\n')

        if is_fastq:
            read_quality = process_readline(in_file.readline(), is_python3)
            line = process_readline(in_file.readline(), is_python3)
            while not line.startswith('@'):
                read_quality += line
                line = process_readline(in_file.readline(), is_python3)
                if not line:
                    if fpos == in_file.tell():
                        break
                    fpos = in_file.tell()
            if len(read_value) != len(read_quality):
                error('The length of sequence and quality lines should be the same! '
                      'Check read %s (SEQ length is %d, QUAL length is %d)' %
                      (read_name, len(read_value), len(read_quality)), log)
            out_file.write('+\n')
            out_file.write(read_quality + '\n')
        return line  # next read name or empty string

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
                        input_file = gzip.open(interlaced_reads, 'r')
                        ungzipped = os.path.splitext(interlaced_reads)[0]
                        out_basename, ext = os.path.splitext(os.path.basename(ungzipped))
                    else:
                        input_file = open(interlaced_reads, 'r')
                        out_basename, ext = os.path.splitext(os.path.basename(interlaced_reads))

                    if interlaced_reads in options_storage.dict_of_prefixes:
                        ext = options_storage.dict_of_prefixes[interlaced_reads]
                    if ext.lower().startswith('.fq') or ext.lower().startswith('.fastq'):
                        is_fastq = True
                        ext = '.fastq'
                    else:
                        is_fastq = False
                        ext = '.fasta'

                    out_left_filename = os.path.join(dst, out_basename + "_1" + ext)
                    out_right_filename = os.path.join(dst, out_basename + "_2" + ext)

                    if not (options_storage.continue_mode and os.path.isfile(out_left_filename) and os.path.isfile(out_right_filename)):
                        options_storage.continue_mode = False
                        log.info("== Splitting " + interlaced_reads + " into left and right reads (in " + dst + " directory)")
                        out_files = [open(out_left_filename, 'w'), open(out_right_filename, 'w')]
                        i = 0
                        next_read_name = write_single_read(input_file, out_files[i], None, is_fastq,
                                                           sys.version.startswith('3.') and was_compressed)
                        while next_read_name:
                            i = (i + 1) % 2
                            next_read_name = write_single_read(input_file, out_files[i], next_read_name, is_fastq,
                                                               sys.version.startswith('3.') and was_compressed)
                        if i == 0:
                            error("The number of reads in file with interlaced reads (" + interlaced_reads + ") should be EVEN!", log)
                        out_files[0].close()
                        out_files[1].close()
                    input_file.close()
                    new_reads_library['left reads'].append(out_left_filename)
                    new_reads_library['right reads'].append(out_right_filename)
                    if interlaced_reads in options_storage.dict_of_prefixes:
                        del options_storage.dict_of_prefixes[interlaced_reads]
                del new_reads_library['interlaced reads']
        new_dataset_data.append(new_reads_library)
    return new_dataset_data


def process_nxmate_reads(dataset_data, dst, log):
    try:
        import lucigen_nxmate
        new_dataset_data = list()
        for reads_library in dataset_data:
            new_reads_library = dict(reads_library)
            if new_reads_library['type'] == 'nxmate':
                raw_left_reads = new_reads_library['left reads']
                raw_right_reads = new_reads_library['right reads']
                new_reads_library['left reads'] = []
                new_reads_library['right reads'] = []
                new_reads_library['single reads'] = []
                for id, left_reads_fpath in enumerate(raw_left_reads):
                    right_reads_fpath = raw_right_reads[id]
                    processed_left_reads_fpath, processed_right_reads_fpath, single_reads_fpath = \
                        lucigen_nxmate.process_reads(left_reads_fpath, right_reads_fpath, dst, log)
                    new_reads_library['left reads'].append(processed_left_reads_fpath)
                    new_reads_library['right reads'].append(processed_right_reads_fpath)
                    new_reads_library['single reads'].append(single_reads_fpath)
                new_reads_library['type'] = 'mate-pairs'
                new_reads_library['orientation'] = 'fr'
            new_dataset_data.append(new_reads_library)
        return new_dataset_data
    except ImportError:
        error("Can't process Lucigen NxMate reads! lucigen_nxmate.py is missing!", log)


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


def read_fasta(filename, gzipped=False):
    res_name = []
    res_seq = []
    first = True
    seq = ''
    if gzipped:
        file_handler = gzip.open(filename)
    else:
        file_handler = open(filename)
    for line in file_handler:
        line = process_readline(line, gzipped and sys.version.startswith('3.'))
        if not line:
            continue
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
    file_handler.close()
    return zip(res_name, res_seq)


def write_fasta(filename, fasta):
    outfile = open(filename, 'w')
    for name, seq in fasta:
        outfile.write(name + '\n')
        for i in range(0, len(seq), 60):
            outfile.write(seq[i : i + 60] + '\n')
    outfile.close()


def break_scaffolds(input_filename, threshold, replace_char="N", gzipped=False):
    new_fasta = []
    modified = False
    for id, (name, seq) in enumerate(read_fasta(input_filename, gzipped)):
        i = 0
        cur_contig_number = 1
        cur_contig_start = 0
        while (i < len(seq)) and (seq.find("N", i) != -1):
            if replace_char != "N":
                modified = True
            start = seq.find("N", i)
            end = start + 1
            while (end != len(seq)) and (seq[end] == "N"):
                end += 1

            i = end + 1
            if (end - start) >= threshold:
                modified = True
                if cur_contig_start != start:
                    new_fasta.append((name.split()[0] + "_" + str(cur_contig_number) + " " + " ".join(name.split()[1:]),
                                      seq[cur_contig_start:start].replace("N", replace_char)))
                    cur_contig_number += 1
                cur_contig_start = end
        if cur_contig_start < len(seq):
            new_fasta.append((name.split()[0] + "_" + str(cur_contig_number) + " " + " ".join(name.split()[1:]),
                              seq[cur_contig_start:].replace("N", replace_char)))
    return modified, new_fasta


def comp(letter):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}[letter.upper()]


def rev_comp(seq):
    return ''.join(itertools.imap(comp, seq[::-1]))


def get_contig_id(s):
    values = s.split("_")
    if len(values) < 2 or (values[0] != ">NODE" and values[0] != "NODE"):
        warning("Contig %s has unknown ID format" % (s))
        return None
    if s.find("'") != -1:
        return (values[1] + "'")
    return values[1]


def remove_fasta_pref(s):
    if s.startswith(">"):
        return s[1:]
    return s


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def is_int(value):
    try:
        int(value)
        return True
    except ValueError:
        return False
