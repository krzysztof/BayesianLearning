import numpy as np
import sympy as sp
import random
import sys
from fractions import Fraction

def _add_bs(b1, b2):
    return b1 + b2
def _sub_bs(b1, b2):
    return b1 - b2
def _mult_b(b, c):
    return b * c
def _div_b(b, c):
    return b / c

def _add_bs_p(b1, b2):
    return b1 * b2
def _sub_bs_p(b1, b2):
    return b1 / b2
def _mult_b_p(b, c):
    return b ** c
def _div_b_p(b, c):
    #print b, c, b**Fraction(1,c)
    return b ** (Fraction(1,c))
    #return b ** (1./c)

SUM_OPERATIONS = (_add_bs, _sub_bs, _mult_b, _div_b, ) 
PRODUCT_OPERATIONS = (_add_bs_p, _sub_bs_p, _mult_b_p, _div_b_p, ) 

def firsti(v,func):
    """Returns index of first element in v meeting func, -1 if otherwise"""
    for i in xrange(len(v)):
        if func(v[i]):
            return i
    return -1

def GaussJordanElimination(vvi_coefficients, vf_absTerms=None, elementary_operations = SUM_OPERATIONS):
    """
        @param vvi_coefficients: coefficients matrix
        @type vvi_coefficients: numpy.matrix (2D)
        @param vf_absTerms: absolute terms
        @type vf_absTerms: numpy.matrix (1D)
        @return: (reduced row echelon form, absolute terms)
        @rtype: (numpy.matrix, numpy.matrix)
    """
    DEBUG = False
    #elementary_operations = SUM_OPERATIONS
    #elementary_operations = PRODUCT_OPERATIONS
    add_bs, sub_bs, mult_b, div_b = elementary_operations
    A = vvi_coefficients.copy()
    b = vf_absTerms.copy() if vf_absTerms != None else None

    pivots = [False]*len(A)

    for col in xrange(A.shape[1]):
        row = firsti([A[:,col][i] if not pivots[i] else 0 for i in xrange(len(pivots))], lambda x: abs(x) > 10e-9)
        if row == -1:
            continue
        pivots[row]=True
        pivot = A[row][col]
        if DEBUG:
            print "Step %d: , pivot: %d" % (col, pivot),

        A[row] /= pivot #Fraction(pivot,1)

        if DEBUG:
            print "row: %s" % (A[row]),
        if b != None:
            #b[row] = b[row]/pivot
            #b[row] = mult_b(b[row], 1/pivot)
            b[row] = div_b(b[row], pivot)
            if DEBUG:
                print "new b: %s" % (b[row])

        for r in xrange(len(A[:,col])):
            if r == row:
                continue
            if abs(A[r][col]) > 10e-7:
                c = A[r][col]
                A[r] = (A[r] - c*A[row]).copy()
                if b != None:
                    #b[r] = sub_bs(b[r], c*b[row])
                    b[r] = sub_bs(b[r], mult_b(b[row],c))

        if DEBUG:
            print "new A: %s" % (A)
            print "new b: %s" % (b)
    return A, b
            

def _get_random_matrix(m,n):
    arr = []
    #print "hi", random.randint(0,30)
    for i in range(m):
        x = random.randint(1,2**n-1)
        while(x in arr):
            x = random.randint(1,2**n-1)
        arr.append(x)
    
    arr = [ [Fraction(y) for y in bin(x)[2:].zfill(n)] for x in arr]
    
    A = np.array(arr)
    return A

def get_random_matrix(m,n, fun = lambda x: True):
    A = _get_random_matrix(m,n)
    while(not fun(A)):
        A = _get_random_matrix(m,n)
    return A

#def print_for_n(N):
#    for i in range(20):
#        A = get_random_matrix(N)
#        while(np.linalg.det(A) != 0):
#            A = get_random_matrix(N)
#        print "\n",A
#A = np.array([[1,1,1,0],[1,0,1,1],[1,1,0,0],[1,0,1,0]])

def augment(matrices):
    return np.column_stack(matrices)

def main():
    #random.seed(1234)
    #A = get_random_matrix(int(sys.argv[1]), fun = lambda x: np.linalg.det(x)!=0)
    #A = np.array([[1,1,0.,1],[1,1,0.,0.],[0.,1,1,0.],[0.,0.,1,0.],[1,0.,0.,0.],[1,0.,0.,1]])
    #A = np.array([[1,1,0,1],[1,1,0,0],[0,1,1,0],[0,0,1,0],[1,0,0,0],[1,0,0,1]])

    A = get_random_matrix(int(sys.argv[1]), int(sys.argv[2]))
    b = np.array([sp.Symbol('b%d'%(i),real=True) for i in range(1, A.shape[0]+1)])
    print augment([A, b])

    print "After"
    A2, b2 = GaussJordanElimination(A, b)
    print augment([A2, b2])
    #print A2
    #print b2.reshape(b2.shape[0],1)

if __name__ == "__main__":
    #print_for_n(int(sys.argv[1]))
    main()
    #pass
