#!/usr/bin/env python
import sys
import getopt
import os
import logging

import haplocontigs_io
import configs_dir

ds_pipeline_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
ds_home = os.path.dirname(os.path.dirname(ds_pipeline_dir))
ds_binary_dir = os.path.join(ds_home, "build/release/bin/")
sp_pipeline_dir = os.path.join(os.path.dirname(ds_pipeline_dir), "spades_pipeline/")

sys.path.insert(0, os.path.abspath(sp_pipeline_dir))
import process_cfg
import support

class DS_Args_List:
	long_options = "allow-gaps weak-align hap=".split()
	short_options = "o:"

class DS_Args:
	allow_gaps = False
	weak_align = False
	haplocontigs_fnames = []
	output_dir = ""
	haplocontigs = ""

	def print_fields(self):
		print "allow_gaps: " + str(self.allow_gaps)
		print "weak_align: " + str(self.weak_align)
		print "files with haplocontigs : " + str(self.haplocontigs_fnames)
		print "haplocontigs file: " + str(self.haplocontigs)
		print "output_dir: " + str(self.output_dir)

def call_method(obj, method_name):
	getattr(obj, method_name)()

def usage():
	print "todo: print usage"

def parse_arguments(argv):
	try:
		options, not_options = getopt.gnu_getopt(argv, DS_Args_List.short_options, DS_Args_List.long_options)
	except getopt.GetoptError as err:
		print str(err) 
        	usage()
        	sys.exit(2)

	ds_args = DS_Args()
	for opt, arg in options:
		if opt == '-o':
			ds_args.output_dir = arg
		elif opt == '--allow-gaps':
			ds_args.allow_gaps = True
		elif opt == '--weak_align':
			ds_args.weak_align = True
		elif opt == '--hap':
			ds_args.haplocontigs_fnames.append(arg)
	ds_args.haplocontigs = os.path.join(ds_args.output_dir, "haplocontigs")
	return ds_args

def check_output_dir(output_dir):
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

def create_log(output_dir):
	log = logging.getLogger('dipspades')
	log.setLevel(logging.DEBUG)
    	console = logging.StreamHandler(sys.stdout)
    	console.setFormatter(logging.Formatter('%(message)s'))
    	console.setLevel(logging.DEBUG)
    	log.addHandler(console)

	log_filename = os.path.join(output_dir, "dipspades.log")
	if os.path.exists(log_filename):
		os.remove(log_filename)
        log_handler = logging.FileHandler(log_filename, mode='a')
   	log.addHandler(log_handler)
	return log

def check_binary(binary_dir, log):
	binary_path = os.path.join(binary_dir, "dipspades")
        if not os.path.isfile(binary_path):
            support.error("DipSPAdes binaries not found: " + binary_path, log)
	return binary_path


def get_dict_of_args(ds_args):
	args_dict = dict()
	args_dict["tails_lie_on_bulges"] = process_cfg.bool_to_str(ds_args.allow_gaps)
	args_dict["align_bulge_sides"] = process_cfg.bool_to_str(ds_args.weak_align)
	args_dict["haplocontigs"] = ds_args.haplocontigs
	args_dict["output_dir"] = ds_args.output_dir
	return args_dict

def prepare_config(config_fname, ds_args, log):
	args_dict = get_dict_of_args(ds_args)
	process_cfg.substitute_params(config_fname, args_dict, log)

def dipspades():
	ds_args = parse_arguments(sys.argv)
	call_method(ds_args, "print_fields")

	check_output_dir(ds_args.output_dir)
	haplocontigs_io.write_haplocontigs_in_file(ds_args.haplocontigs, ds_args.haplocontigs_fnames)
	print "file with haplocontigs was written\n"

	log = create_log(ds_args.output_dir)
	print "log was created"

	config_fname = configs_dir.prepare_configs(os.path.join(ds_home, "configs/dipspades"), ds_args)
	prepare_config(config_fname, ds_args, log)
	print "configs were prepared"

	binary_path = check_binary(ds_binary_dir, log)
	command = [binary_path, config_fname]
	support.sys_call(command, log)

dipspades()
