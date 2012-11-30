import unittest

from bayesiandataset import BayesianDataSet, ChildNode

class BayesianDataSetTest(unittest.TestCase):
	def setUp(self):
		self.bn0 = BayesianDataSet('testData/CancerMAX1.txt') # THIS IS NOT 20k test case
		self.bn0.addChildNode(4,'No', [0,1,2,3], ['False', 'False', 'False', 'Medium'])
		self.child_node0 = self.bn0.children[4]

	def test_encoding(self):
		#print self.bn0.domain
		#print self.bn0._str_to_num
		#print self.child_node0

		#domain is:
		#[ 0 :[0,1,2], 1:[0, 1], 2:[0, 1], 3:[0, 1, 2]]
		# parent char states are:
		# [2,1,1,1]
		# string domain:
		##(('OnePack', 'TwoPacks', 'False'), ('True', 'False'), ('True', 'False'), ('Bad', 'Medium', 'Good'), ('Benign', 'Malignant', 'No'))

		e0 = self.bn0._decode_equation([1,0,0,0,1,0], self.child_node0)
		print [self.bn0.domain[i][e0[i]] for i in range(len(e0))]
		pairs = (
			((1,0,0,1), (1,0,1,1,0,0)),
			((2,1,0,0), (0,0,0,1,0,1)),
			((2,0,0,0), (0,0,1,1,0,1)),
			((2,1,1,1), (0,0,0,0,0,0)),
			((0,0,0,0), (0,1,1,1,0,1)),
			# single params
			((1,1,1,1), (1,0,0,0,0,0)),
			((0,1,1,1), (0,1,0,0,0,0)),
			((2,0,1,1), (0,0,1,0,0,0)),
			((2,1,0,1), (0,0,0,1,0,0)),
			((2,1,1,2), (0,0,0,0,1,0)),
			((2,1,1,0), (0,0,0,0,0,1)),
		)
		for orig, enc in pairs:
			enc0 = self.bn0._encode_equation(orig, self.child_node0)
			self.assertEqual(tuple(enc0), enc)

			dec0 = self.bn0._decode_equation(enc, self.child_node0)
			self.assertEqual(tuple(dec0), orig)

def main():
	unittest.main()

if __name__ == '__main__':
	main()