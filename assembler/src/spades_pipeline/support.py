#!/usr/bin/env python

import os
import stat
import sys

# Based on http://stackoverflow.com/a/616686/92396
class Tee(object):
    def __init__(self, name, mode, console=True):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        self.console = console
        sys.stdout = self
        sys.stderr = self

    def free(self):
        sys.stdout = self.stdout
        sys.stderr = self.stderr
        self.file.close()

    def write(self, data):
        self.file.write(data)
        if self.console:
            self.stdout.write(data)
        self.flush()

    def flush(self):
        self.file.flush()
        self.stdout.flush()


def verify(expr, message):
    if (not (expr)):
        print ("Assertion failed. Message: " + message)
        exit(1)


def error(err_str, prefix="== Error == "):
    print >> sys.stderr, "\n\n" + prefix + " " + err_str + "\n\n"
    exit(1)


def warning(warn_str, prefix="== Warning == "):
    print("\n\n" + prefix + " " + warn_str + "\n\n")


#TODO: error log -> log
#TODO: os.sytem gives error -> stop

def sys_call(cmd, cwd=None):
    import shlex
    import time
    import sys
    import subprocess

    cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd)

    while not proc.poll():
        if proc.returncode is not None:
            break
        sys.stdout.write(proc.stdout.readline())
        time.sleep(0)

    for line in proc.stdout.readlines():
        print line,
    proc.communicate()

    if proc.returncode:
        error('system call for: "%s" finished abnormally, err code: %d' % (cmd, proc.returncode))


def save_data_to_file(data, file):
    output = open(file, 'wb')
    output.write(data.read())
    output.close()
    os.chmod(file, stat.S_IWRITE | stat.S_IREAD | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
