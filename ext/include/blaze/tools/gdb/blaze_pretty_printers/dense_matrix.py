""" Pretty printers for dense matrix classes
"""

from .util import stripType, templateParams, strToBool, dataPtrFromAlignedArray


class DenseMatrixPrinter:
	"""Basic class for dense matrix printers"""

	def children(self):
		"""Enumerate child items"""
		
		yield ('rows', self.rows)
		yield ('columns', self.columns)
		yield ('data', self.data)

		# Enumerate elements
		for i in range(self.rows):
			for j in range(self.columns):
				ind = j * self.spacing + i if self.columnMajor else i * self.spacing + j
				yield ('[{0},{1}]'.format(i, j), self.data[ind])
		

	def storageOrder(self):
		return 'columnMajor' if self.columnMajor else 'rowMajor'


	# def display_hint(self):
	# 	return 'array'


class StaticMatrixPrinter(DenseMatrixPrinter):
	"""Print blaze::StaticMatrix"""

	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		template_params = templateParams(self.type)

		self.rows = int(template_params[1])
		self.columns = int(template_params[2])
		self.columnMajor = strToBool(template_params[3])
		self.spacing = val['MM'] if self.columnMajor else val['NN']		
		self.elementType = self.type.template_argument(0)
		self.val = val
		self.data = dataPtrFromAlignedArray(self.val['v_'])

		
	def to_string(self):
		return 'StaticMatrix<{0}, {1}, {2}, {3}>'.format(
			self.elementType, self.rows, self.columns, self.storageOrder())


class HybridMatrixPrinter(DenseMatrixPrinter):
	"""Print blaze::HybridMatrix"""

	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		template_params = templateParams(self.type)

		self.rows = val['m_']
		self.columns = val['n_']
		self.maxRows = int(template_params[1])
		self.maxColumns = int(template_params[2])
		self.columnMajor = strToBool(template_params[3])
		self.spacing = val['MM'] if self.columnMajor else val['NN']
		self.elementType = self.type.template_argument(0)		
		self.val = val
		self.data = dataPtrFromAlignedArray(self.val['v_'])

		
	def to_string(self):
		return 'HybridMatrix<{0}, {1}, {2}, {3}>'.format(
			self.elementType, self.maxRows, self.maxColumns, self.storageOrder())


class DynamicMatrixPrinter(DenseMatrixPrinter):
	"""Print blaze::DynamicMatrix"""

	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		template_params = templateParams(self.type)
		
		self.rows = val['m_']
		self.columns = val['n_']
		self.columnMajor = strToBool(template_params[1])
		self.spacing = val['mm_'] if self.columnMajor else val['nn_']		
		self.elementType = self.type.template_argument(0)		
		self.val = val
		self.data = self.val['v_']
			
	
	def to_string(self):
		return 'DynamicMatrix<{0}, {1}>'.format(self.elementType, self.storageOrder())


class CustomMatrixPrinter(DenseMatrixPrinter):
	"""Print blaze::CustomMatrix"""

	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		template_params = templateParams(self.type)
		
		self.rows = val['m_']
		self.columns = val['n_']

		self.aligned = strToBool(template_params[1])
		self.padded = strToBool(template_params[2])
		self.columnMajor = strToBool(template_params[3])
		self.spacing = val['mm_'] if self.columnMajor else val['nn_']		
		self.elementType = self.type.template_argument(0)		
		self.val = val
		self.data = self.val['v_']


	def alignment(self):
		return 'aligned' if self.aligned else 'unaligned'


	def padding(self):
		return 'padded' if self.padded else 'unpadded'
			
	
	def to_string(self):
		return 'CustomMatrix<{0}, {1}, {2}, {3}>'.format(
			self.elementType, self.alignment(), self.padding(), self.storageOrder())
