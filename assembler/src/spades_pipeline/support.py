#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import stat
import sys


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


#TODO: error log -> log
#TODO: os.sytem gives error -> stop

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


def save_to_yaml(data, filename):
    INDENT = '  '
    yaml = open(filename, 'w')
    cur_indent = 0
    yaml.write(cur_indent * INDENT + '[')
    cur_indent += 1
    yaml.write('\n' + cur_indent * INDENT + '{')
    cur_indent += 1
    first = True
    for key, value in data.iteritems():
        if not first:
            yaml.write(',')
        yaml.write('\n' + cur_indent * INDENT + key + ': ')
        if isinstance(value, list):
            yaml.write('[')
            cur_indent += 1
            first = True
            for v in value:
                if not first:
                    yaml.write(',')
                yaml.write('\n' + cur_indent * INDENT + '"' + str(v) + '"')
                first = False
            cur_indent -= 1
            yaml.write('\n' + cur_indent * INDENT + ']')
        else:
            yaml.write('"' + str(value) + '"')
        first = False
    cur_indent -= 1
    yaml.write('\n' + cur_indent * INDENT + '}')
    cur_indent -= 1
    yaml.write('\n' + cur_indent * INDENT + ']')


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