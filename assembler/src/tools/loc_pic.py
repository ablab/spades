#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import string
import re
import subprocess
import datetime

init_dir = "before_simplification/pos_loc"
dirs = ["tip_clipping_", "err_con_removal_", "bulge_removal_"]
add_dirs = ["final_tip_clipping", "final_bulge_removal", "final_err_con_removal", "final_simplified"]
try:
	g_pos = subprocess.check_output('ls -1 ' + init_dir, shell=True)
except:
	print "Not found:", directory
	exit(0)
g_pos = g_pos.strip().split('\n')

output_root = "interest_pos/"
os.mkdir(output_root)

for pos in g_pos:
	os.mkdir(output_root + pos)

cnt = 0;
for i in range(0, 10):
	for dir_pref in dirs:
		if (os.path.exists(dir_pref + str(i))):
			print "found dir " + dir_pref + str(i) 
			for pos in g_pos:
				file_from = dir_pref + str(i) + "/pos_loc/" + pos + "/kmer1_.dot"
				print "looking for file " + file_from
				if (os.path.exists(file_from)):
					print "found file " + file_from
					file_to = output_root + pos + "/" + str(cnt) + "_" + dir_pref + str(i) + ".dot"
					print "copying to " + file_to
					os.system('cp ' + file_from + " " + file_to) 
			cnt = cnt + 1	
		else:
			print "Finished"
			exit(0)

