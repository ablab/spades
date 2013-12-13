#!/usr/bin/env python

import os
import sys

class ConfigField:
	name = ""
	value = ""
	comment = ";"

	def __init__(self, new_name, new_value, new_comment = ""):
		self.name = new_name
		self.value = new_value
		self.comment = new_comment

class Subconfig:
	field_list = list()
	name = ""

	def __init__(self, new_name):
		self.name = new_name

class DSConfig:
	subconfigs = list()
	fields = list()


