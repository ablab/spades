""" Pretty printers for sparse matrix classes
"""

from .util import stripType, templateParams, strToBool


class CompressedMatrixPrinter:
	"""Print blaze::StaticMatrix"""

	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		template_params = templateParams(self.type)

		self.rows = val['m_']
		self.columns = val['n_']
		self.begin = val['begin_']
		self.end = val['end_']
		self.columnMajor = strToBool(template_params[1])
		self.elementType = self.type.template_argument(0)
		

	def children(self):
		"""Enumerate child items"""	
		
		yield ('rows', self.rows)
		yield ('columns', self.columns)
		
		# Enumerate non-zero elements and their indices
		for i in range(self.columns if self.columnMajor else self.rows):
			begin = self.begin[i]
			end = self.end[i]
			
			it = begin
			while it < end:
				j = it.dereference()['index_']
				row = j if self.columnMajor else i
				column = i if self.columnMajor else j
				yield ('[{0}, {1}]'.format(row, column), it.dereference()['value_'])
				it += 1


	def storageOrder(self):
		return 'columnMajor' if self.columnMajor else 'rowMajor'

		
	def to_string(self):
		return "CompressedMatrix<{0}, {1}>".format(self.elementType, self.storageOrder())
