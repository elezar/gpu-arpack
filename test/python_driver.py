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
	func = lib['dense_seev']
	func.argtypes = [ct.c_int, ct.POINTER(ct.c_float), 
					 ct.c_int, ct.c_int, 
					 ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
					 ct.POINTER(ct.c_double), ct.POINTER(ct.c_int)]

	info_fp = np.zeros(10, dtype=ct.c_double)
	info_int = np.zeros(10, dtype=ct.c_int)

	N = 100
	LDA = 100

	NEV = 10

	A = np.diag(np.arange(1,N+1),k=0).astype(ct.c_float)

	eigenvalues = np.zeros(NEV, dtype=np.complex64)
	eigenvectors = np.zeros((N,NEV), dtype=ct.c_float)


	func(N, as_pointer(A, ct.c_float), 
		LDA, NEV, 
		as_pointer(eigenvalues, ct.c_float), as_pointer(eigenvectors, ct.c_float), 
		as_pointer(info_fp, ct.c_double), as_pointer(info_int, ct.c_int))


	print eigenvalues

	max_eig = N

	print '%5s:' % 'index', '%15s' % ('calculated'), '%15s' % ('reference'), '%15s' % 'relative error'
	for i in range(len(eigenvalues)):
		ref = max_eig-i
		err = np.abs(eigenvalues[i] - ref)/np.abs(ref)
		print '%5d:' % (i+1), '%15s' % (eigenvalues[i]), '%15s' % (ref), '%15.3e' % err



if __name__ == '__main__':
	main()
