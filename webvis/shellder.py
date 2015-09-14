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
        prepare_pipe(pipe_in)
        prepare_pipe(pipe_out)
        self.pout = open(pipe_out, "r+")
        self.end_out = end_out
        self.pin = open(pipe_in, "w+")
        pin = open(pipe_in, "r+")
        pout = open(pipe_out, "w+")
        pushd = os.getcwd()
        #Launch the reader
        self.queue = Queue.Queue()
        self.reader = None
        #Launch the process
        os.chdir(dir)
        self.proc = subprocess.Popen(["./run", "rv"], stdin=pin, stdout=pout)
        os.chdir(pushd)

    def get_line(self, timeout=None):
        try:
            return self.queue.get(True, timeout)
        except Queue.Empty:
            print "No output"
            return None

    #Reads the whole output and returns as a list of strings
    def get_output(self, timeout=None):
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
        res = []
        while True:
            str = self.get_line(timeout)
            if str is None:
                break
            res.append(str)
            if str == self.end_out:
                break
        return res

    #Sends a command, like "help"
    def send(self, command):
        self.pin.write(command + "\n")
        self.pin.flush()
        return self

    def close(self):
        self.proc.kill()
        self.pout.close()
        self.pin.close()
