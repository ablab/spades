#!/usr/bin/env python

import os
import sys
import shutil

def check_config_directory(output_base):
	config_dir = os.path.join(output_base, "configs")
	#if not os.path.exists(config_dir):
	#	os.makedirs(config_dir)
	return config_dir

# src_config_dir - path of dipspades configs
def copy_configs(src_config_dir, dst_config_dir):
	if os.path.exists(dst_config_dir):
		shutil.rmtree(dst_config_dir)
	shutil.copytree(src_config_dir, dst_config_dir)

def prepare_configs(src_config_dir, ds_args):
	config_dir = check_config_directory(ds_args.output_dir)
	print "Config directory " + config_dir + "\n"
	copy_configs(src_config_dir, config_dir)
	print "Configs were copied"
	config_fname = os.path.join(config_dir, "config.info")
	return os.path.abspath(config_fname)
