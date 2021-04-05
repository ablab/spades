""" Pretty printers for dense vector classes
"""

from .util import stripType, templateParams, strToBool, dataPtrFromAlignedArray


class DenseVectorPrinter:
	"""Basic class for dense vector printers"""

	def children(self):
		"""Enumerate child items"""

		yield ('size', self.size)
		yield ('data', self.data)

		# Enumerate elements
		for i in range(self.size):
			yield ('[{0}]'.format(i), self.data[i])


	def transposeFlag(self):
		return 'rowVector' if self.rowVector else 'columnVector'


	# def display_hint(self):
	# 	return 'array'


class StaticVectorPrinter(DenseVectorPrinter):
	"Print blaze::StaticVector"

	def __init__(self, val):
		"Extract all the necessary information"

		self.type = stripType(val.type)
		template_params = templateParams(self.type)

		self.size = int(template_params[1])		
		self.rowVector = strToBool(template_params[2])		
		self.elementType = self.type.template_argument(0)	
		self.val = val
		self.data = dataPtrFromAlignedArray(self.val['v_'])

		
	def to_string(self):
		return "StaticVector<{0}, {1}, {2}>".format(
			self.elementType, self.size, self.transposeFlag())


class HybridVectorPrinter(DenseVectorPrinter):
	"Print blaze::HybridVector"

	def __init__(self, val):
		"Extract all the necessary information"

		self.type = stripType(val.type)
		template_params = templateParams(self.type)

		self.size = val['size_']
		self.maxSize = int(template_params[1])
		self.rowVector = strToBool(template_params[2])		
		self.elementType = self.type.template_argument(0)		
		self.val = val
		self.data = dataPtrFromAlignedArray(self.val['v_'])

		
	def to_string(self):
		return "HybridVector<{0}, {1}, {2}>".format(
			self.elementType, self.maxSize, self.transposeFlag())


class DynamicVectorPrinter(DenseVectorPrinter):
	"Print blaze::DynamicVector"

	def __init__(self, val):
		"Extract all the necessary information"

		self.type = stripType(val.type)
		template_params = templateParams(self.type)
		
		self.size = val['size_']		
		self.rowVector = strToBool(template_params[1])		
		self.elementType = self.type.template_argument(0)
		self.val = val
		self.data = self.val['v_']
			
	
	def to_string(self):
		return "DynamicVector<{0}, {1}>".format(self.elementType, self.transposeFlag())


class CustomVectorPrinter(DenseVectorPrinter):
	"Print blaze::CustomVector"

	def __init__(self, val):
		"Extract all the necessary information"

		self.type = stripType(val.type)
		template_params = templateParams(self.type)
		
		self.size = val['size_']
		self.aligned = strToBool(template_params[1])
		self.padded = strToBool(template_params[2])
		self.rowVector = strToBool(template_params[3])		
		self.elementType = self.type.template_argument(0)		
		self.val = val
		self.data = self.val['v_']


	def alignment(self):
		return 'aligned' if self.aligned else 'unaligned'


	def padding(self):
		return 'padded' if self.padded else 'unpadded'
			
	
	def to_string(self):
		return "CustomVector<{0}, {1}, {2}, {3}>".format(
			self.elementType, self.alignment(), self.padding(), self.transposeFlag())
