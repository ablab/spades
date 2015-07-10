import os
import subprocess
from time import sleep

#Making in-out pipes
def prepare_pipe(path, mode):
    if not os.path.exists(path):
        os.mkfifo(path)
    return open(path, mode)

pipe_in = prepare_pipe("/tmp/vis_in", "r+")   #for sending commands
pipe_out = prepare_pipe("/tmp/vis_out", "w+") #for receiving output

#Hold pipe
#TODO: check if necessary
pipe_dontclose = open("/tmp/vis_out", "r+")

#Starting the visualizer
print os.getcwd()
os.chdir("../../../") #TODO: softcode
print os.getcwd()
proc = subprocess.Popen(["./run", "rv"], stdin=pipe_in, stdout=pipe_out)

#Wait eternally
while True:
    sleep(10000)
