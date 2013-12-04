#!/usr/bin/env python

import os
import sys

import config_structure

def is_subconfig(string):
	print "impement me!"
	return True

def is_field(string):
	print "impement me!"
	return True

def get_subconfig_name(name_string):
	splits = name_string.split()
	if len(splits) == 0:
		return ""
	s_index = 0
	for s in splits:
		index = find(s, "{")
		if index != -1:
			if len(s) == 1:
				if s_index == 0:
					return ""
				return splits[s_index - 1]
			return s[:index]
		s_index = s_index + 1
	return ""

def is_subconfig_end(string):
	splits = string.split()
	if length(splits) == 1:
		if splits[0] == "}":
			return True
	return False

def parse_field(string):
	splits = string.split()
	comment = ""
	comment_starts = False
	name = splits[0]
	value = splits[1]
	for s in splits:
		if find(s, ";") != -1:
			comment = s
			comment_starts = True
		elif comment_starts:
			comment = comment + s
	field = config_structure.ConfigField(name, value, comment)
	return field

def parse_subconfig(name_string, reader):
	subconfig = config_structure.Subconfig(get_subconfig_name(name_string))
	for line in reader:
		if is_subconfig_end(line):
			break
		elif is_field(line):
			new_field = parse_field(line)
			subconfig.field_list.append(new_field)
	return subconfig

def read_config_from_file(filename):
	ds_config = config_structure.DSConfig()
	reader = open(filename, "r")
	for line in reader:
		if is_subconfig(line):
			new_subconfig = parse_subconfig(line, reader)
			ds_config.subconfigs.append(new_subconfig)
		#elif is_field(line):
		#	new_field = parse__field(line)
		#	ds_config.fields.append(new_field)
	return ds_config
