import os
import subprocess
import sys
from time import sleep

pipe_in = open(sys.argv[1], "r+")  #for sending commands
pipe_out = open(sys.argv[2], "w+") #for receiving output

#Hold pipe
#TODO: check if necessary
pipe_dontclose = open(sys.argv[2], "r+")

#Starting the visualizer
os.chdir(sys.argv[3]) #TODO: softcode
proc = subprocess.Popen(["./run", "rv"], stdin=pipe_in, stdout=pipe_out)

#Wait eternally
while True:
    sleep(10000)
