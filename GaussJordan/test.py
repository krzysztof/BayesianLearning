import unittest
from fractions import Fraction
from __init__ import augment, GaussJordanElimination, PRODUCT_OPERATIONS, SUM_OPERATIONS
import sympy as sp
import numpy as np

def cmp_matrices(A, B):
    assert A.shape == B.shape, "Matrices of different size."
    item_cnt = reduce(lambda x,y: x*y, A.shape)
    return all(A.reshape(item_cnt)[i] == B.reshape(item_cnt)[i] for i in range(item_cnt))

class TestSolving(unittest.TestCase):
    def setUp(self):
        self.A1_raw = [[-3, 1, 1, 1, 1],
            [-2, 1, 1, 1, 0],
            [-2, 1, 0, 1, 1],
            [-2, 0, 1, 1, 1],
            [-2, 1, 1, 0, 1]]
        self.A1 = np.array(self.A1_raw)*Fraction(1,1)
        self.b1 = np.array([sp.Symbol('b%d'%(i),real=True) for i in range(self.A1.shape[0])])

    def testSum1(self):
        RREF, b2 =  GaussJordanElimination(self.A1, self.b1, elementary_operations = SUM_OPERATIONS)
        b = self.b1
        b2_correct = np.array([-3*b[0] + b[1] + b[2] + b[3] + b[4],
                    -2*b[0] + b[1] + b[2] + b[4],
                    -2*b[0] + b[1] + b[3] + b[4],
                    -2*b[0] + b[1] + b[2] + b[3],
                    -2*b[0] + b[2] + b[3] + b[4]])
        assert cmp_matrices(RREF, np.eye(5,dtype=object))
        assert cmp_matrices(b2, b2_correct)

if __name__ == "__main__":
    unittest.main()
