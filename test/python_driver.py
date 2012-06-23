__author__ = "Evan Lezar"
__date__ = "17 June 2012"

"""A quick test for the ARPACK eigensolver.
"""

import numpy as np
import ctypes as ct

def as_pointer(A, ctype_type = None):
    return A.ctypes.data_as(ct.POINTER(ctype_type))

def main():
	lib = ct.CDLL('../lib/libsingle_dense_cpu.so')
	func = lib['dense_sgev']
	func.argtypes = [ct.c_int, ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), 
					 ct.c_int, ct.c_int, ct.c_float, ct.POINTER(ct.c_float), 
					 ct.POINTER(ct.c_float),
					 ct.POINTER(ct.c_double), ct.POINTER(ct.c_int)]

	info_fp = np.zeros(10, dtype=ct.c_double)
	info_int = np.zeros(10, dtype=ct.c_int)

	N = 100
	LDA = 100

	NEV = 10

	S = np.diag(np.arange(1,N+1),k=0).astype(ct.c_float)
	T = np.eye(N, N, dtype=ct.c_float)

	shift = np.zeros(1, dtype=ct.c_float)

	eigenvalues = np.zeros(NEV, dtype=np.complex64)
	eigenvectors = np.zeros((N,NEV), dtype=ct.c_float)


	func(N, as_pointer(S, ct.c_float), as_pointer(T, ct.c_float), 
		LDA, NEV, shift, 
		as_pointer(eigenvalues, ct.c_float), as_pointer(eigenvectors, ct.c_float), 
		as_pointer(info_fp, ct.c_double), as_pointer(info_int, ct.c_int))


	print eigenvalues



if __name__ == '__main__':
	main()
