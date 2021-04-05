""" Pretty printers for sparse vector classes
"""

from .util import stripType, templateParams, strToBool


class CompressedVectorPrinter:
	"Print blaze::CompressedVector"

	def __init__(self, val):
		"""Extract all the necessary information"""

		self.type = stripType(val.type)
		template_params = templateParams(self.type)

		self.rowVector = strToBool(template_params[1])		
		self.elementType = self.type.template_argument(0)	
		self.val = val
		self.size = val['size_']
		self.begin = val['begin_']
		self.end = val['end_']


	def children(self):
		"""Enumerate child items"""

		yield ('size', self.size)

		# Enumerate non-zero elements and their indices
		cur = self.begin
		while cur < self.end:
			yield ('[{0}]'.format(cur.dereference()['index_']), cur.dereference()['value_'])
			cur += 1


	def transposeFlag(self):
		return 'rowVector' if self.rowVector else 'columnVector'

		
	def to_string(self):
		return "CompressedVector<{0}, {1}>".format(self.elementType, self.transposeFlag())
