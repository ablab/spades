
############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from os import path
import os
import Queue
import subprocess
import threading

#Making in-out pipes
def prepare_pipe(path):
    if not os.path.exists(path):
        os.mkfifo(path)

class Shellder:
    #Initialized for sending commands to the vis, reading output, or both
    def __init__(self, pipe_in, pipe_out, dir="", end_out="[end]\n"):
        #Open two-sided named pipes
        self.pipe_in = pipe_in
        self.pipe_out = pipe_out
        self.end_out = end_out
        #Prepare the reader
        self.queue = Queue.Queue()
        self.proc = None
        self.reader = None
        self.dir = dir
        #Open IO pipes
        prepare_pipe(pipe_in)
        prepare_pipe(pipe_out)
        self.cout = open(pipe_out, "w+")
        self.pout = open(pipe_out, "r+")
        self.cin = open(pipe_in, "r+")
        self.pin = open(pipe_in, "w+")
        self.launch_if_needed()

    def launch_if_needed(self):
        #Launch online_vis process if it hasn't started yet or was aborted
        if self.proc is None or self.proc.poll() is not None:
            #Launch the process
            #pin = open(self.pipe_in, "r+")
            #pout = open(self.pipe_out, "w+")
            self.proc = subprocess.Popen(["./run", "rv"], cwd=self.dir,
                 stdin=self.cin, stdout=self.cout, stderr=self.cout)

        #Launch the reader thread
        if self.reader is None or not self.reader.is_alive():
            def reader():
                #TODO: refactor doublecode
                while True:
                    str = self.pout.readline()
                    self.queue.put(str)
                    if str == self.end_out:
                        break
            self.reader = threading.Thread(target=reader)
            self.reader.start()

    def get_line(self, timeout=None):
        try:
            return self.queue.get(True, timeout)
        except Queue.Empty:
            return None

    #Reads the whole output and returns as a list of strings
    def get_output(self, timeout=None):
        res = []
        complete = False
        try:
            self.launch_if_needed()
            while True:
                str = self.get_line(timeout)
                if str is None:
                    break
                if str == self.end_out:
                    complete = True
                    break
                if str.startswith("Exception"):
                    raise OSError(str)
                res.append(str)
        except IOError as e:
            res.append("Cannot communicate with online_vis: {}\n".format(str(e)))
            complete = True
        if timeout is None:
            return res
        return (res, complete)

    #Sends a command, like "help"
    def send(self, command):
        self.pin.write(command + "\n")
        self.pin.flush()
        return self

    def close(self):
        self.proc.kill()
        self.pout.close()
        self.pin.close()

    def pid(self):
        return self.proc.pid
