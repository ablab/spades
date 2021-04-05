""" Pretty printer for blaze::Vector class
"""

from .util import stripType, templateParams, strToBool, dataPtrFromAlignedArray


class VectorPrinter:
	"""Printer for blaze::Vector"""
	
	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		self.derivedType = self.type.template_argument(0)
		self.derived = val.cast(self.derivedType)


	def children(self):
		"""Enumerate elements"""

		yield (str(self.derivedType), self.derived)
