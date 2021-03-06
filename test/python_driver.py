__author__ = "Evan Lezar"
__date__ = "17 June 2012"

"""A quick test for the ARPACK eigensolver.
"""

import numpy as np
import ctypes as ct

def as_pointer(A, ctype_type = None):
    return A.ctypes.data_as(ct.POINTER(ctype_type))

def main():
	# Load the shared library.
	lib = ct.CDLL('../lib/libsingle_dense_cpu.so')

	# Select the function to call.
	# This has the prototype
	# dense_seev(int, float*, int, int, float*, float*, double*, int*)
	func = lib['dense_seev']
	func.argtypes = [ct.c_int, ct.POINTER(ct.c_float), 
					 ct.c_int, ct.c_int, 
					 ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),
					 ct.POINTER(ct.c_double), ct.POINTER(ct.c_int)]

	# Allocate storage for timing and info.
	info_fp = np.zeros(10, dtype=ct.c_double)
	info_int = np.zeros(10, dtype=ct.c_int)

	# Setup the eigenvalue problem.
	N = 100
	LDA = N
	A = np.diag(np.arange(1,N+1),k=0).astype(ct.c_float)

	# The number of eigenvalues and storage for the eigenvalues and eigenvectors.
	NEV = 10
	eigenvalues = np.zeros(NEV, dtype=np.complex64)
	eigenvectors = np.zeros((N,NEV), dtype=ct.c_float)

	# call the loaded C function.
	func(N, as_pointer(A, ct.c_float), 
		LDA, NEV, 
		as_pointer(eigenvalues, ct.c_float), as_pointer(eigenvectors, ct.c_float), 
		as_pointer(info_fp, ct.c_double), as_pointer(info_int, ct.c_int))

	# Output the results.
	max_eig = N
	print 'Results:'
	print '%5s:' % 'index', '%15s' % ('calculated'), '%15s' % ('reference'), '%15s' % 'relative error'
	for i in range(len(eigenvalues)):
		ref = max_eig-i
		err = np.abs(eigenvalues[i] - ref)/np.abs(ref)
		print '%5d:' % (i+1), '%15s' % (eigenvalues[i]), '%15s' % (ref), '%15.3e' % err

	print 'Timing:'
	print info_fp
	print info_int

	
if __name__ == '__main__':
	main()
