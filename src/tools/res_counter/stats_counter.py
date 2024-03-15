#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import shutil
import re
import time
import commands

timeout = 10

def get_memory_stats(pid):
    pstree = os.path.join(os.path.abspath(sys.path[0]), "./pstree")
    output = commands.getoutput(pstree + " " + str(pid))
    lines = output.splitlines()
    
    mem_used = 0
    for line in lines:
        if line.find('%') == -1:
            continue
        cur_pid = int(line[line.find('%') + 2:].split()[0])    #     |     \-+-% 27271 gurevich FindErrors
        output = commands.getoutput("ps aux | grep ' " + str(cur_pid) + " '").splitlines() ##
        for finded in output:
            if (finded.split()[1] == str(cur_pid)):
                mem_used += int(finded.split()[5])
                break	

    return mem_used


if len(sys.argv) < 2:
	print 'ERROR!'
	print 'Usage:', sys.argv[0], 'pid'
	sys.exit(0)

pid = sys.argv[1]

#memory stats
max_mem_used = 0

#time stats
elapsed = ""

#try:
while (os.path.exists("/proc/" + str(pid))): 
  cur_mem_used = get_memory_stats(pid)    
  if (max_mem_used  < cur_mem_used):
    max_mem_used  = cur_mem_used
  #print output
  output = commands.getoutput("ps -p " + str(pid) + " -o etime")
  if (len(output.splitlines()) > 1):
    elapsed = output.splitlines()[1].strip()
  #print output
  time.sleep(timeout)
#except:
#   pass

###
print >> sys.stderr, "Total stats:"
print >> sys.stderr, "time elapsed: ", elapsed
print >> sys.stderr, "max mem_used: %db %dM %.1fG" % (max_mem_used * 1024, max_mem_used // 1024, max_mem_used // 1024 / 1024.0)
