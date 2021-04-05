""" Pretty printer for blaze::Matrix class
"""

from .util import stripType, templateParams, strToBool, dataPtrFromAlignedArray


class MatrixPrinter:
	"""Printer for blaze::Matrix"""
	
	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		self.derivedType = self.type.template_argument(0)
		self.derived = val.cast(self.derivedType)


	def children(self):
		"""Enumerate elements"""

		yield (str(self.derivedType), self.derived)
