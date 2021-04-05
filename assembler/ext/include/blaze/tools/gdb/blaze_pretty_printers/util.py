""" Utility functions
"""

import gdb
import re


def stripType(type):
	"""Returns type with references, qualifiers, and typedefs removed"""
	
	if type.code == gdb.TYPE_CODE_REF:
		type = type.target()

	return type.unqualified().strip_typedefs()


def templateParams(type):
	"""Return template parameters as a list of strings

	The gdb extension does not support value template arguments -- need to extract them by hand.
	"""

	m = re.match('[\w:]+<(.*)>', type.tag)
	if m:
		template_params = re.findall('[\w:]+(?:<.*>)?', m[1])
		template_params = [x.replace(" ", "") for x in template_params]
		return template_params
	else:
		return None


def strToBool(s):
	"""Convert 'true' or 'false' string to a bool value
	"""

	if s == 'true':
		return True
	elif s == 'false':
		return False
	else:
		raise ValueError('Boolean value must be "true" or "false"')


def dataPtrFromAlignedArray(aligned_array):
	"""Get data pointer from an AlignedArray object
	"""
	element_type = aligned_array.type.template_argument(0)
	return aligned_array['v_'].cast(element_type.pointer())