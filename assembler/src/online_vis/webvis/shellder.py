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
        def reader(q):
            while True:
                q.put(self.pout.readline())
        self.reader = threading.Thread(target=reader, args=(self.queue,))
        self.reader.start()
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
        res = []
        str = self.get_line(timeout)
        while str is not None:
            res.append(str)
            if str == self.end_out:
                break
            str = self.get_line(timeout)
        return res

    #Sends a command, like "help"
    def send(self, command):
        self.pin.write(command + "\n")
        self.pin.flush()
        return self

    def close(self):
        #TODO: fix close error
        self.reader.join(0.1)
        self.pin.close()
        self.pout.close()
        self.proc.kill()
