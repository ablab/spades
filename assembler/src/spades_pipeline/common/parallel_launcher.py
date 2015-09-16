############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import multiprocessing
import os
import shlex
import subprocess
import sys
import signal
import time
import traceback
import logging

__author__ = 'anton'
import subprocess

class FunctionTask:
    def __init__(self, f):
        self.f = f

    def __call__(self, barcode):
        try:
            self.f(barcode)
        except:
            traceback.print_tb(sys.exc_info()[2])
            return 1
        return 0

class PseudoLambda:
    def __init__(self):
        pass
    def __call__(self, task):
        task.run()

def GetHandlers(output_file_pattern, err_file_pattern, bid):
    if output_file_pattern == "":
        output_file_pattern = "/dev/null"
    output = open(output_file_pattern.format(bid), "a")
    if err_file_pattern == "":
        return (output, subprocess.STDOUT)
    else:
        return (output, open(err_file_pattern.format(bid), "a"))

class ExternalCallTask:
    def __init__(self, output_pattern = "", err_pattern = "", log_name = None):
        self.output_pattern = output_pattern
        self.err_pattern = err_pattern
        self.log_name = log_name

    def __call__(self, data):
        bid, command = data
        output, err = GetHandlers(self.output_pattern, self.err_pattern, bid)
        if self.log_name is not None:
            logging.getLogger(self.log_name).info("Starting: " + str(bid))
        import shlex
        return_code = subprocess.call(shlex.split(command), stdout = output, stderr = err)
        if return_code == 0:
            logging.getLogger(self.log_name).info("Successfully finished: " + str(bid))
        else:
            logging.getLogger(self.log_name).info("Failed to finish: " + str(bid))
        return return_code



def run_in_parallel(task, material, threads):
    result = call_in_parallel(task, material, threads)
    errors = len(material) - result.count(0)
    return errors

def call_in_parallel(task, material, threads):
    pool = multiprocessing.Pool(threads)
    result = pool.map_async(task, material).get(1000000000)
#    result = pool.map(call, commands)
    return result
