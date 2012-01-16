#!/usr/bin/python

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

for pos in g_pos
	os.mkdir(output_root + g_pos)

cnt = 1;
for i in range(1, 10)
	for dir_pref in dirs
		if (os.path.exists(dir_pref + str(i)))
			for pos in g_pos
				file_from = dir_pref + str(i) + "/" + pos
				if (os.path.exists(file_from))
					file_to = output_root + pos + "/" + cnt + "_" + dir_pref + str(i)
					os.system('cp ' + file_from + " " + file_to) 
			cnt++		
		else
			print "Finished"
			exit(0)

