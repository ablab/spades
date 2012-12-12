#!/usr/bin/env python

import os
import stat
import sys

def verify(expr, log, message):
    if (not (expr)):
        log.info ("Assertion failed. Message: " + message)
        sys.exit(1)


def error(err_str, log, prefix="== Error == "):
    log.info("\n\n" + prefix + " " + err_str + "\n\n")
    sys.exit(1)


def warning(warn_str, log, prefix="== Warning == "):
    log.info("\n\n" + prefix + " " + warn_str + "\n\n")


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


def save_data_to_file(data, file):
    output = open(file, 'wb')
    output.write(data.read())
    output.close()
    os.chmod(file, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
