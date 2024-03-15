""" Pretty printers for Blaze classes.

Usage:

Create a ~/.gdbinit file with the following content:
	python
	import sys
	sys.path.insert(0, '/path/to/blaze/tools/gdb')
	import blaze_pretty_printers
	end

For more information on writing gdb pretty printers,
please refer to https://sourceware.org/gdb/onlinedocs/gdb/Pretty-Printing-API.html
"""

import gdb
import re
import itertools

from .dense_vector import StaticVectorPrinter, HybridVectorPrinter, DynamicVectorPrinter, CustomVectorPrinter
from .dense_matrix import StaticMatrixPrinter, HybridMatrixPrinter, DynamicMatrixPrinter, CustomMatrixPrinter
from .sparse_vector import CompressedVectorPrinter
from .sparse_matrix import CompressedMatrixPrinter
from .vector import VectorPrinter
from .matrix import MatrixPrinter
from .util import stripType


pretty_printers_dict = {
	re.compile('^blaze::StaticMatrix<.*>$'): lambda val: StaticMatrixPrinter(val),
	re.compile('^blaze::HybridMatrix<.*>$'): lambda val: HybridMatrixPrinter(val),
	re.compile('^blaze::DynamicMatrix<.*>$'): lambda val: DynamicMatrixPrinter(val),
	re.compile('^blaze::CustomMatrix<.*>$'): lambda val: CustomMatrixPrinter(val),
	re.compile('^blaze::StaticVector<.*>$'): lambda val: StaticVectorPrinter(val),
	re.compile('^blaze::HybridVector<.*>$'): lambda val: HybridVectorPrinter(val),
	re.compile('^blaze::DynamicVector<.*>$'): lambda val: DynamicVectorPrinter(val),
	re.compile('^blaze::CustomVector<.*>$'): lambda val: CustomVectorPrinter(val),
	re.compile('^blaze::CompressedVector<.*>$'): lambda val: CompressedVectorPrinter(val),
	re.compile('^blaze::CompressedMatrix<.*>$'): lambda val: CompressedMatrixPrinter(val),
	re.compile('^blaze::Matrix<.*>$'): lambda val: MatrixPrinter(val),
	re.compile('^blaze::Vector<.*>$'): lambda val: VectorPrinter(val),
}


def lookup_function(val):
	"""Look-up and return a pretty-printer that can print val.
	"""
	
	type = stripType(val.type)
	
	typename = type.tag
	if typename == None:
		return None
	
	for function in pretty_printers_dict:
		if function.search(typename):
			return pretty_printers_dict[function](val)
	
	return None


# Register the lookup function with GDB
gdb.pretty_printers.append(lookup_function)
