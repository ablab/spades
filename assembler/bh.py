#!/usr/bin/env python

import os
import sys
import datetime
import shutil
import subprocess

sys.path.append("src/tools/spades_pipeline/")

import support

from support import redirected_stream

from process_cfg import *

sys.path.append("src/tools")
from dataset import hammer

def prepare_config(filename, cfg, build_path):

	subst_dict = dict()

	subst_dict["input_file_0"]            = cfg.left_reads
	subst_dict["input_file_1"]            = cfg.right_reads
	subst_dict["input_working_dir"]  	  = build_path
	subst_dict["general_max_iterations"]  = cfg.iteration_count 
	subst_dict["input_qvoffset"]          = cfg.qvoffset 

	substitute_params(filename, subst_dict)

def main():
	#Functions
	def build_folder(cfg):
		import datetime
		suffix = datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
		return cfg.working_dir + r"/" + cfg.dataset_name + r"/" + suffix + r"/"
	
	def build_bh(cfg):
		support.sys_call(cfg.output_to_console, cfg.log_filename, "make -C ./build/release/hammer")

	def copy(dest):	
		os.makedirs(dest, 0755)
		shutil.copy("src/hammer/config.inp", dest)
		shutil.copy("build/release/hammer/main", dest)
	
	def verify(expr, message):
		if (not (expr)):
			print "Assertion failed. Message: " + message
			exit(0)

	def check_files(files, str_it_count):
		msg = "Failure! Check output files."
		verify(len(files) == 4, msg + "1")
		verify(str_it_count + ".reads.left.corrected" in files, msg + "2")
		verify(str_it_count + ".reads.right.corrected" in files, msg + "3")
		verify(str_it_count + ".reads.left.unpaired" in files, msg + "4")
		verify(str_it_count + ".reads.right.unpaired" in files, msg + "5")

	def determine_read_files(folder, str_it_count):
		check_files(subprocess.check_output('ls -1 ' + folder + r"/" + str_it_count + "*.reads.* | xargs -n1 basename"\
		, shell=True).strip().split('\n'), str_it_count)
		
		answer = dict()	
		answer["first"] = folder + str_it_count + ".reads.left.corrected"	
		answer["second"] = folder + str_it_count + ".reads.right.corrected"
		answer["single_first"] = folder + str_it_count + ".reads.left.unpaired"	
		answer["single_second"] = folder + str_it_count + ".reads.right.unpaired"
		return answer

	def generate_dataset(cfg, tmp_dir):	
		str_it_count = str(cfg.iteration_count - 1)
		if (cfg.iteration_count <= 10):
			str_it_count = "0" + str_it_count
		dataset_cfg = determine_read_files(tmp_dir + r"/", str_it_count)
		dataset_cfg["RL"] = "100"
		dataset_cfg["single_cell"] = "true"
		return cfg.dataset_name + "\n" + hammer(dataset_cfg, cfg.output_dir, True)
		
	
	#Preparations
	cfg = load_config_from_file(config_file_name())
	build_path = build_folder(cfg)
	print("Path " + build_path)
	os.makedirs(build_path, 0755)

	log_filename = build_path + "bayes_hammer.log"
	cfg.__dict__["log_filename"] = log_filename

	print("\n== Log can be found here: " + cfg.log_filename + "\n")

	log_file = open(log_filename, "w")

	old_stdout = sys.stdout
	old_stderr = sys.stderr

	if cfg.output_to_console:
		sys.stderr = redirected_stream(log_file, sys.stderr)
		sys.stdout = redirected_stream(log_file, sys.stdout)
	else:
		sys.stderr = redirected_stream(log_file, None)
		sys.stdout = redirected_stream(log_file, None)

	#Process started
	print("\n== Compilation started. Don't start another instance of BayesHammer before application is started ==\n")  
	build_bh(cfg)
	print("\n== Compilation finished successfully. ==\n")

	copy(build_path + "exec")
	prepare_config(build_path + "exec/config.inp", cfg, build_path)
	print("\n== Starting BayesHammer. You can start new instance now. ==\n") 
    
	support.sys_call(cfg.output_to_console, cfg.log_filename, build_path + "exec/main " + build_path + "exec/config.inp")

	sys.stdout = old_stdout
	sys.stderr = old_stderr
	
	print("\n== BayesHammer run finished. ==\n") 
	print("\n== Generating dataset " + cfg.dataset_name + ". ==\n") 
	
	#os.system("touch " + build_path + "02.reads.left.corrected "  + build_path + "02.reads.left.unpaired "  + build_path + "02.reads.right.corrected " + build_path + "02.reads.right.unpaired")

	dataset_str = generate_dataset(cfg, build_path)
	
	datasets_file = open("configs/debruijn/datasets.info", "a")	
	datasets_file.write("\n" + dataset_str)

	print(dataset_str)

	print("\n== Script finished. ==\n") 

if __name__ == '__main__':
	main()

