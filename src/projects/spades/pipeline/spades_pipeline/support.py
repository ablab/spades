#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import gzip
import logging
import os
import shutil
import stat
import sys
from platform import uname


log = logging.getLogger("spades")

# constants to print and detect warnings and errors in logs
SPADES_PY_ERROR_MESSAGE = "== Error == "
SPADES_PY_WARN_MESSAGE = "== Warning == "
SPADES_ERROR_MESSAGE = " ERROR "
SPADES_WARN_MESSAGE = " WARN "
# for correct warnings detection in case of continue_mode
continue_logfile_offset = None
# for removing tmp_dir even if error occurs

only_old_style_options = True
old_style_single_reads = False


def error(err_str, logger_instance=None, prefix=SPADES_PY_ERROR_MESSAGE, exit_code=-1):
    binary_name = "SPAdes"

    if logger_instance:
        logger_instance.error("\n\n%s %s" % (prefix, err_str))
        log_warnings(logger_instance, with_error=True)
        logger_instance.info("\nIn case you have troubles running %s, you can report an issue on our GitHub repository "
                             "github.com/ablab/spades" % binary_name)
        logger_instance.info(
            "Please provide us with params.txt and %s.log files from the output directory." % binary_name.lower())
    else:
        sys.stderr.write("\n\n%s %s\n\n" % (prefix, err_str))
        sys.stderr.write(
            "\nIn case you have troubles running %s, you can report an issue on our GitHub repository "
            "github.com/ablab/spades\n" % binary_name)
        sys.stderr.write(
            "Please provide us with params.txt and %s.log files from the output directory.\n" % binary_name.lower())
        sys.stderr.flush()
    if current_tmp_dir and os.path.isdir(current_tmp_dir):
        shutil.rmtree(current_tmp_dir)
    sys.exit(exit_code)


def warning(warn_str, logger_instance=None, prefix=SPADES_PY_WARN_MESSAGE):
    if logger_instance:
        logger_instance.warning("\n\n%s %s\n\n" % (prefix, warn_str))
    else:
        sys.stdout.write("\n\n%s %s\n\n\n" % (prefix, warn_str))
        sys.stdout.flush()


def wsl_check():
    def in_wsl():
        return 'microsoft' in uname()[2].lower()

    if in_wsl():
        return ("1. WSL is an unsupported platform\n"
                "2. If SPAdes crashes, then you might want to compile it from sources\n"
                "3. If nothing works, run on real Linux")
    return ""


def get_error_hints(exit_code):
    if exit_code == -11:
        return wsl_check()


def sys_error(cmd, logger_instance, exit_code):
    hints_str = get_error_hints(exit_code)
    err_msg = "system call for: \"%s\" finished abnormally, OS return value: %d\n%s" % (cmd, exit_code, hints_str)
    error(err_msg, logger_instance, exit_code=exit_code)


# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    def is_exe(file_path):
        return os.path.isfile(file_path) and os.access(file_path, os.X_OK)

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
                    avail_mem = int(line[len(avail_mem_header):].split()[0])  # in kB
                    avail_mem /= 1024 * 1024  # in GB
                    return avail_mem
        except ValueError:
            return None
        except IOError:
            return None
    return None


# based on http://stackoverflow.com/questions/196345/how-to-check-if-a-string-in-python-is-in-ascii
def is_ascii_string(line):
    try:
        line.encode("ascii")
    except UnicodeDecodeError:  # python2
        return False
    except UnicodeEncodeError:  # python3
        return False
    else:
        return True


def process_readline(line, utf8=True):
    if utf8:
        return str(line, "utf-8").rstrip()
    return line.rstrip()


def process_spaces(string):
    if " " in string:
        string = '"' + string + '"'
    return string


def sys_call(cmd, logger_instance=None, cwd=None):
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    output = ""
    while not proc.poll():
        line = process_readline(proc.stdout.readline())
        if line:
            if logger_instance:
                logger_instance.info(line)
            else:
                output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = process_readline(line)
        if line:
            if logger_instance:
                logger_instance.info(line)
            else:
                output += line + "\n"

    if proc.returncode:
        sys_error(cmd, logger_instance, proc.returncode)
    return output


def universal_sys_call(cmd, logger_instance, out_filename=None, err_filename=None, cwd=None):
    # Runs cmd and redirects stdout to out_filename (if specified), stderr to err_filename (if specified), or to log otherwise
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

    if logger_instance and (not out_filename or not err_filename):
        while not proc.poll():
            if not out_filename:
                line = process_readline(proc.stdout.readline())
                if line:
                    logger_instance.info(line)
            if not err_filename:
                line = process_readline(proc.stderr.readline())
                if line:
                    logger_instance.info(line)
            if proc.returncode is not None:
                break

        if not out_filename:
            for line in proc.stdout.readlines():
                if line != "":
                    logger_instance.info(process_readline(line))
        if not err_filename:
            for line in proc.stderr.readlines():
                if line != "":
                    logger_instance.info(process_readline(line))
    else:
        proc.wait()

    if out_filename:
        stdout.close()
    if err_filename:
        stderr.close()
    if proc.returncode:
        sys_error(cmd, logger_instance, proc.returncode)


def save_data_to_file(data, file):
    with open(file, "wb") as output:
        output.write(data.read())

    os.chmod(file, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def get_important_messages_from_log(log_filename, warnings=True):
    def already_saved(list_to_check, suffix_str):  # for excluding duplicates (--continue-from may cause them)
        for item in list_to_check:
            if item.endswith(suffix_str):
                return True
        return False

    if warnings:
        spades_py_message = SPADES_PY_WARN_MESSAGE
        spades_message = SPADES_WARN_MESSAGE
    else:  # errors
        spades_py_message = SPADES_PY_ERROR_MESSAGE
        spades_message = SPADES_ERROR_MESSAGE

    # for capturing correct warnings in case of continue_mode
    if continue_logfile_offset:
        with open(log_filename) as continued_log:
            continued_log.seek(continue_logfile_offset)
            continued_stage_phrase = continued_log.readline()
            while not continued_stage_phrase.strip():
                continued_stage_phrase = continued_log.readline()
            lines_to_check = continued_log.readlines()

        all_lines = open(log_filename).readlines()
        failed_stage_index = all_lines.index(continued_stage_phrase)
        lines_to_check = all_lines[:failed_stage_index] + lines_to_check
    else:
        lines_to_check = open(log_filename).readlines()

    spades_py_msgs = []
    spades_msgs = []
    important_message_summary_prefix = " * "
    for line in lines_to_check:
        if line.startswith(important_message_summary_prefix):
            continue
        if spades_py_message in line:
            suffix = line[line.find(spades_py_message) + len(spades_py_message):].strip()
            line = line.replace(spades_py_message, "").strip()
            if not already_saved(spades_py_msgs, suffix):
                spades_py_msgs.append(important_message_summary_prefix + line)
        elif spades_message in line:
            suffix = line[line.find(spades_message) + len(spades_message):].strip()
            line = line.strip()
            if not already_saved(spades_msgs, suffix):
                spades_msgs.append(important_message_summary_prefix + line)
    return spades_py_msgs, spades_msgs


def get_logger_filename(logger_instance):
    log_file = None
    for h in logger_instance.__dict__["handlers"]:
        if h.__class__.__name__ == "FileHandler":
            log_file = h.baseFilename
    return log_file


def log_warnings(logger_instance, with_error=False):
    log_file = get_logger_filename(logger_instance)
    if not log_file:
        return False
    for h in logger_instance.__dict__["handlers"]:
        h.flush()
    spades_py_warns, spades_warns = get_important_messages_from_log(log_file, warnings=True)
    if spades_py_warns or spades_warns:
        if with_error:
            logger_instance.warning("\n======= SPAdes pipeline finished abnormally and WITH WARNINGS!")
        else:
            logger_instance.info("\n======= SPAdes pipeline finished WITH WARNINGS!")
        warnings_filename = os.path.join(os.path.dirname(log_file), "warnings.log")
        warnings_handler = logging.FileHandler(warnings_filename, mode='w')
        logger_instance.addHandler(warnings_handler)
        # log.info("===== Warnings occurred during SPAdes run =====")
        logger_instance.info("")
        if spades_py_warns:
            logger_instance.info("=== Pipeline warnings:")
            for line in spades_py_warns:
                logger_instance.info(line)
        if spades_warns:
            logger_instance.info("=== Error correction and assembling warnings:")
            for line in spades_warns:
                logger_instance.info(line)
        logger_instance.info("======= Warnings saved to " + warnings_filename)
        logger_instance.removeHandler(warnings_handler)
        if with_error:
            spades_py_errors, spades_errors = get_important_messages_from_log(log_file, warnings=False)
            logger_instance.info("")
            logger_instance.info("=== ERRORs:")
            for line in (spades_errors + spades_py_errors):
                logger_instance.info(line)
        return True
    return False


# START for processing YAML files


def pretty_print_reads(dataset_data, indent="    "):
    read_types = ["left reads", "right reads", "interlaced reads", "single reads", "merged reads"]
    for lib_id, reads_library in enumerate(dataset_data):
        log.info(indent + "Library number: %d, library type: %s" % (lib_id + 1, reads_library["type"]))
        if "orientation" in reads_library:
            log.info("%s  orientation: %s" % (indent, reads_library["orientation"]))
        for reads_type in read_types:
            if reads_type not in reads_library:
                value = "not specified"
            else:
                value = str(reads_library[reads_type])
            log.info("%s  %s: %s" % (indent, reads_type, value))

# END: for processing YAML files


def read_fasta(filename, gzipped=False):
    res_name = []
    res_seq = []
    first = True
    seq = ""
    if gzipped:
        file_handler = gzip.open(filename)
    else:
        file_handler = open(filename)
    for line in file_handler:
        line = process_readline(line, gzipped)
        if not line:
            continue
        if line[0] == '>':
            res_name.append(line.strip())
            if not first:
                res_seq.append(seq)
            else:
                first = False
            seq = ""
        else:
            seq += line.strip()
    res_seq.append(seq)
    file_handler.close()
    return zip(res_name, res_seq)


def write_fasta(filename, fasta):
    with open(filename, 'w') as outfile:
        for name, seq in fasta:
            outfile.write(name + '\n')
            for i in range(0, len(seq), 60):
                outfile.write(seq[i: i + 60] + '\n')


def break_scaffolds(input_filename, threshold, replace_char="N", gzipped=False):
    new_fasta = []
    modified = False
    for lib_id, (name, seq) in enumerate(read_fasta(input_filename, gzipped)):
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
    return "".join(map(comp, seq[::-1]))


def get_contig_id(s):
    values = s.split("_")
    if len(values) < 2 or (values[0] != ">NODE" and values[0] != "NODE"):
        warning("contig %s has unknown ID format" % s)
        return None
    if "'" in s:
        return values[1] + "'"
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


def add_user_write_permission_recursive(path):
    s = os.stat(path)
    os.chmod(path, s.st_mode | stat.S_IWUSR | stat.S_IRUSR)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            d = os.path.join(root, d)
            s = os.stat(d)
            os.chmod(d, s.st_mode | stat.S_IWUSR | stat.S_IRUSR)
        for f in files:
            f = os.path.join(root, f)
            s = os.stat(f)
            os.chmod(f, s.st_mode | stat.S_IWUSR | stat.S_IRUSR)


# shutil.copyfile does not copy any metadata (time and permission), so one
# cannot expect preserve_mode = False and preserve_times = True to work.
def copy_tree(src, dst, preserve_times=True, preserve_mode=True):
    if not preserve_mode:
        copy_fn = shutil.copyfile
    else:
        copy_fn = shutil.copy2

    if os.path.exists(dst):
        shutil.rmtree(dst)

    # shutil.copytree preserves the timestamp, so we must update it afterwards.
    shutil.copytree(src, dst, copy_function=copy_fn)
    if not preserve_mode:
        add_user_write_permission_recursive(dst)

    if not preserve_times:
        for dirpath, _, filenames in os.walk(dst):
            os.utime(dirpath, None)
            for file in filenames:
                os.utime(os.path.join(dirpath, file), None)
