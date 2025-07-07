#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import glob
import itertools
import logging
import math
import os
import re
import shutil
import sys
import tempfile
import traceback
from os.path import abspath, expanduser, join

from .common import SeqIO
from .options_storage import OptionStorage
from .support import error, is_ascii_string, warning, only_old_style_options, old_style_single_reads

options_storage = OptionStorage()


log = logging.getLogger("spades")

def check_python_version():
    if sys.version_info < options_storage.MINIMAL_PYTHON_VERSION:
        error("\nPython version %s is not supported!\n"
              "Minimal supported version is %s" %
              (sys.version.split()[0], ".".join(list(map(str, options_storage.MINIMAL_PYTHON_VERSION)))))
        return False
    return True


def get_spades_binaries_info_message():
    return "You can obtain SPAdes binaries in one of two ways:" + \
           "\n1. Download them from https://github.com/ablab/spades/" + \
           "\n2. Build source code with ./spades_compile.sh script"


def check_binaries(binary_dir):
    for binary in ["spades-hammer", "spades-ionhammer", "spades-core", "spades-bwa"]:
        binary_path = os.path.join(binary_dir, binary)
        if not os.path.isfile(binary_path):
            error("SPAdes binaries not found: %s\n%s" % (binary_path, get_spades_binaries_info_message()), log)


def check_file_existence(input_filename, message="", logger_instance=None):
    filename = abspath(expanduser(input_filename))
    check_path_is_ascii(filename, message)
    if not os.path.isfile(filename):
        error("file not found: %s (%s)" % (filename, message), logger_instance=logger_instance)
    options_storage.dict_of_rel2abs[input_filename] = filename
    return filename


def get_read_file_type(input_filename, logger_instance=None):
    if input_filename in options_storage.dict_of_prefixes:
        ext = options_storage.dict_of_prefixes[input_filename]
        file_type = SeqIO.get_read_file_type("filename" + ext)
    else:
        file_type = SeqIO.get_read_file_type(input_filename)

    if not file_type:
        error("incorrect extension of reads file: %s" % input_filename, logger_instance)
    return file_type


def check_file_not_empty(input_filename, message="", logger_instance=None):
    filename = abspath(expanduser(input_filename))
    file_type = get_read_file_type(input_filename, logger_instance)
    if file_type == 'bam' or file_type == 'sra':
        return

    try:
        reads_iterator = SeqIO.parse(SeqIO.Open(filename, "r"), file_type)
        if next(reads_iterator, None) is None:
            error("file is empty: %s (%s)" % (filename, message), logger_instance=logger_instance)
    except Exception as inst:
        error(inst.args[0].format(FILE=filename) + "\n\n" +
              traceback.format_exc().format(FILE=filename), logger_instance=logger_instance)


def check_dir_existence(input_dirname, message="", logger_instance=None):
    dirname = abspath(expanduser(input_dirname))
    check_path_is_ascii(dirname, message)
    if not os.path.isdir(dirname):
        error("directory not found: %s (%s)" % (dirname, message), logger_instance=logger_instance)
    options_storage.dict_of_rel2abs[input_dirname] = dirname
    return dirname


def check_path_is_ascii(path, message=""):
    if not is_ascii_string(path):
        error("path contains non-ASCII characters: %s (%s)" % (path, message))


def recreate_dir(dirname):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)


def check_files_duplication(filenames):
    for filename in filenames:
        if filenames.count(filename) != 1:
            error("file %s was specified at least twice" % filename, log)


def check_reads_file_format(filename, message, only_assembler, iontorrent, library_type):
    if filename in options_storage.dict_of_prefixes:
        ext = options_storage.dict_of_prefixes[filename]
    else:
        ext = os.path.splitext(filename)[1]
        if ext.lower() == ".gz":
            pre_ext = os.path.splitext(filename[:-len(ext)])[1]
            if (pre_ext + ext).lower() in options_storage.ALLOWED_READS_EXTENSIONS:
                ext = pre_ext + ext
            else:  # allows ".fastq.1.gz" like extensions
                pre_pre_ext = os.path.splitext(filename[:-len(pre_ext + ext)])[1]
                ext = pre_pre_ext + ext
    if ext.lower() not in options_storage.ALLOWED_READS_EXTENSIONS:
        error("file with reads has unsupported format (only %s are supported): %s (%s)" %
              (", ".join(options_storage.ALLOWED_READS_EXTENSIONS), filename, message), log)

    if not iontorrent and ext.lower() in options_storage.IONTORRENT_ONLY_ALLOWED_READS_EXTENSIONS:
        error(", ".join(options_storage.IONTORRENT_ONLY_ALLOWED_READS_EXTENSIONS) +
              " formats supported only for iontorrent mode: %s (%s)" % (filename, message), log)

    if (not only_assembler and ext.lower() not in options_storage.BH_ALLOWED_READS_EXTENSIONS and
            library_type not in options_storage.LONG_READS_TYPES):
        error("to run read error correction, reads should be in FASTQ format (%s are supported): %s (%s)" %
              (", ".join(options_storage.BH_ALLOWED_READS_EXTENSIONS), filename, message), log)

    if library_type.endswith("contigs") and ext.lower() not in options_storage.CONTIGS_ALLOWED_READS_EXTENSIONS:
        error("file with %s should be in FASTA format  (%s are supported): %s (%s)" %
              (library_type, ", ".join(options_storage.CONTIGS_ALLOWED_READS_EXTENSIONS), filename, message), log)

    if library_type.endswith("graph") and ext.lower() not in options_storage.GRAPH_ALLOWED_READS_EXTENSIONS:
        error("file with %s should be in GFA format  (%s are supported): %s (%s)" %
              (library_type, ", ".join(options_storage.GRAPH_ALLOWED_READS_EXTENSIONS), filename, message), log)


def get_latest_dir(pattern):
    def atoi(text):
        if text.isdigit():
            return int(text)
        return text

    def natural_keys(text):
        return [atoi(c) for c in re.split("(\\d+)", text)]

    latest_dir = None
    for dir_to_test in sorted(glob.glob(pattern), key=natural_keys, reverse=True):
        if os.path.isdir(dir_to_test):
            latest_dir = dir_to_test
            break
    return latest_dir


def get_tmp_dir(prefix="", base_dir=None):
    global current_tmp_dir
    if not base_dir:
        base_dir = options_storage.args.tmp_dir
    if not os.path.isdir(base_dir):
        os.makedirs(base_dir)
    current_tmp_dir = tempfile.mkdtemp(dir=base_dir, prefix=prefix)
    return current_tmp_dir


def get_short_reads_type(option):
    for short_reads_type in options_storage.SHORT_READS_TYPES.keys():
        if option.startswith("--" + short_reads_type):
            # additional check to except collisions with LONG_READS_TYPES, e.g. --s<#> and --sanger
            if option[len("--" + short_reads_type):len("--" + short_reads_type) + 1].isdigit():
                return short_reads_type
    return None


def get_long_reads_type(option):
    for long_reads_type in options_storage.LONG_READS_TYPES:
        if option.startswith("--") and option in ("--" + long_reads_type):
            return long_reads_type
    return None


def get_graph_type(option):
    for graph_reads_type in options_storage.GRAPH_READS_TYPES:
        if option.startswith("--") and option in ("--" + graph_reads_type):
            return graph_reads_type
    return None


def is_single_read_type(option):
    return option.startswith("--s") and option[3:].isdigit()


def get_lib_type_and_number(option):
    # defaults for simple -1, -2, -s, --12 options
    lib_type = "pe"
    lib_number = 1

    if get_short_reads_type(option):
        lib_type = get_short_reads_type(option)
        lib_number = int(re.search(r'\d+', option).group())
    elif get_long_reads_type(option):
        lib_type = get_long_reads_type(option)
    elif get_graph_type(option):
        lib_type = get_graph_type(option)

    return lib_type, lib_number


def get_data_type(option):
    if option.endswith("-12"):
        data_type = "interlaced reads"
    elif option.endswith("-1"):
        data_type = "left reads"
    elif option.endswith("-2"):
        data_type = "right reads"
    elif option.endswith("-s") or is_single_read_type(option) or get_long_reads_type(option) or get_graph_type(option):
        data_type = "single reads"
    elif option.endswith("-m") or option.endswith("-merged"):
        data_type = "merged reads"
    else:  # -rf, -ff, -fr
        data_type = "orientation"
    return data_type


def get_option_prefix(data):
    prefix = None
    if ':' in data and ('.' + data[:data.find(':')]) in options_storage.ALLOWED_READS_EXTENSIONS:
        prefix = data[:data.find(':')]
        data = data[data.find(':') + 1:]
    return data, prefix


def add_to_dataset(option, data, dataset_data):
    lib_type, lib_number = get_lib_type_and_number(option)
    record_id = "%s_%d" % (lib_type, lib_number)
    data_type = get_data_type(option)
    if data_type == "orientation":
        data = option[-2:]

    if record_id not in dataset_data:  # setting default values for a new record
        dataset_data[record_id] = {}
        dataset_data[record_id]["number"] = lib_number
        if lib_type in options_storage.SHORT_READS_TYPES:
            dataset_data[record_id]["type"] = options_storage.SHORT_READS_TYPES[lib_type]
        else:
            dataset_data[record_id]["type"] = lib_type
    if data_type.endswith("reads"):
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
    for reads_library in dataset_data.values() if isinstance(dataset_data, dict) else dataset_data:
        if not reads_library:
            continue
        has_reads = False
        has_paired_reads = False
        for key in reads_library.keys():
            if key.endswith("reads"):
                has_reads = True
            if key in ["interlaced reads", "merged reads", "left reads", "right reads"]:
                has_paired_reads = True
                break
        if not has_reads:
            continue
        if not has_paired_reads and reads_library["type"] == "paired-end":
            reads_library["type"] = "single"
            if "orientation" in reads_library:
                del reads_library["orientation"]
        if "orientation" not in reads_library:
            if reads_library["type"] == "paired-end" or reads_library["type"] == "hq-mate-pairs":
                reads_library["orientation"] = "fr"
            elif reads_library["type"] == "mate-pairs":
                reads_library["orientation"] = "rf"
        corrected_dataset_data.append(reads_library)
    return corrected_dataset_data


def relative2abs_paths(dataset_data, dirname):
    dirname = abspath(expanduser(dirname))
    abs_paths_dataset_data = []
    for reads_library in dataset_data:
        for key, value in reads_library.items():
            if key.endswith("reads"):
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


def get_reads_length(dataset_data, ignored_types,
                     used_types=options_storage.READS_TYPES_USED_IN_CONSTRUCTION,
                     num_checked=10 ** 4, diff_len_allowable=25):
    max_reads_lenghts = [get_max_reads_length(reads_file, num_checked) for reads_file in
                         get_reads_files(dataset_data, ignored_types, used_types)]

    avg_len = sum(max_reads_lenghts) / len(max_reads_lenghts)
    for max_len in max_reads_lenghts:
        if math.fabs(max_len - avg_len) > diff_len_allowable:
            warning("read lengths differ more than allowable. Length: %f. Avg. length: %f." % (max_len, avg_len), log)
    reads_length = min(max_reads_lenghts)
    log.info("\nReads length: %d\n" % reads_length)
    return reads_length


def get_primary_max_reads_length(dataset_data, ignored_types, used_types, num_checked=10 ** 4):
    max_reads_lenghts = [get_max_reads_length(reads_file, num_checked) for reads_file in
                         get_reads_files(dataset_data, ignored_types, used_types)]

    reads_length = max(max_reads_lenghts)
    log.info("\nReads length: %d\n" % reads_length)
    return reads_length


def get_reads_files(dataset_data, ignored_types, used_types=None):
    for reads_library in dataset_data:
        if (used_types is not None) and reads_library["type"] not in used_types:
            continue
        for key, value in reads_library.items():
            if key in ignored_types:
                log.info("Files with %s were ignored." % key)
                continue
            elif key.endswith("reads"):
                for reads_file in value:
                    yield reads_file


def get_max_reads_length(reads_file, num_checked):
    file_type = get_read_file_type(reads_file, log)
    max_reads_length = 0
    if file_type == 'sra':
        max_reads_length = 100
    else:
        try:
            max_reads_length = max(
                [len(rec) for rec in itertools.islice(SeqIO.parse(SeqIO.Open(reads_file, "r"), file_type), num_checked)])
        except Exception as inst:
            error(inst.args[0].format(FILE=reads_file) + "\n\n" +
                  traceback.format_exc().format(FILE=reads_file), logger_instance=log)
        else:
            log.info("%s: max reads length: %s" % (reads_file, str(max_reads_length)))
    return max_reads_length


def check_dataset_reads(dataset_data, only_assembler, iontorrent):
    all_files = []
    for lib_id, reads_library in enumerate(dataset_data):
        left_number = 0
        right_number = 0
        if "number" not in reads_library:
            reads_library["number"] = lib_id + 1

        for key, value in reads_library.items():
            if key.endswith("reads"):
                for reads_file in value:
                    check_file_existence(reads_file,
                                         "%s, library number: %d, library type: %s" %
                                         (key, reads_library["number"], reads_library["type"]), log)
                    check_reads_file_format(reads_file, "%s, library number: %d, library type: %s" %
                                            (key, reads_library["number"], reads_library["type"]), only_assembler, iontorrent,
                                            reads_library["type"])
                    if reads_library["type"] in options_storage.READS_TYPES_USED_IN_CONSTRUCTION:
                        check_file_not_empty(reads_file,
                                             "%s, library number: %d, library type: %s" %
                                             (key, reads_library["number"], reads_library["type"]), log)

                    all_files.append(reads_file)
                if key == "left reads":
                    left_number = len(value)
                elif key == "right reads":
                    right_number = len(value)

        if left_number != right_number:
            error("the number of files with left paired reads is not equal to the number of files "
                  "with right paired reads (library number: %d, library type: %s)!" %
                  (lib_id + 1, reads_library["type"]), log)
    if not len(all_files):
        error("you should specify at least one file with reads!", log)
    check_files_duplication(all_files)


def check_single_reads_in_options():
    if not only_old_style_options and old_style_single_reads:
        warning("it is recommended to specify single reads with --pe<#>-s, --mp<#>-s, --hqmp<#>-s, "
                "or --s<#> option instead of -s!", log)


def get_lib_ids_by_type(dataset_data, types):
    if type(types) is not list:
        types = [types]
    lib_ids = []
    for lib_id, reads_library in enumerate(dataset_data):
        if reads_library["type"] in types:
            lib_ids.append(lib_id)
    return lib_ids


def get_libs_by_type(dataset_data, types):
    ids = get_lib_ids_by_type(dataset_data, types)
    result = []
    for lib_id in ids:
        result.append(dataset_data[lib_id])
    return result


def rm_libs_by_type(dataset_data, types):
    ids = get_lib_ids_by_type(dataset_data, types)
    for lib_id in sorted(ids, reverse=True):
        del dataset_data[lib_id]
    return dataset_data


def dataset_is_empty(dataset_data):
    for reads_library in dataset_data:
        if reads_library:
            return False
    return True


def dataset_has_gzipped_reads(dataset_data):
    for reads_library in dataset_data:
        for key in reads_library:
            if key.endswith("reads"):
                for reads_file in reads_library[key]:
                    if reads_file.endswith(".gz"):
                        return True
    return False


def dataset_has_interlaced_reads(dataset_data):
    for reads_library in dataset_data:
        if "interlaced reads" in reads_library:
            return True
    return False


def dataset_has_additional_contigs(dataset_data):
    for reads_library in dataset_data:
        if reads_library["type"].endswith("contigs"):
            return True
    return False
