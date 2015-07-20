from os import path
import os
import signal
import subprocess

#Making in-out pipes
def prepare_pipe(path):
    if not os.path.exists(path):
        os.mkfifo(path)

class Shellder:
    #Initialized for sending commands to the vis, reading output, or both
    def __init__(self, pipe_in, pipe_out, dir="", end_out="[end]\n"):
        prepare_pipe(pipe_in)
        prepare_pipe(pipe_out)
        self.pout = open(pipe_out, "r+")
        self.end_out = end_out
        self.pin = open(pipe_in, "w+")
        pin = open(pipe_in, "r+")
        pout = open(pipe_out, "w+")
        os.chdir(dir)
        self.proc = subprocess.Popen(["./run", "rv"], stdin=pin, stdout=pout)

    def get_line(self):
        str = self.pout.readline()
        if str != self.end_out:
            return str
        else:
            return None
        
    #Reads the whole output and returns as a list of strings
    def get_output(self):
        res = []
        str = self.get_line()
        while str is not None:
            res.append(str)
            str = self.get_line()
        return res

    #Sends a command, like "help"
    def send(self, command):
        self.pin.write(command + "\n")
        self.pin.flush()
        return self

    def close(self):
        self.pin.close()
        self.pout.close()
        self.proc.kill()

