""" In this script we will see how to compute matrix - matrix products
using numpy and scipy.linalg.BLAS 

From the output we see the python methods are similar.
"""

from numpy import matrix, array, dot, zeros
from numpy.random import random
from scipy.linalg.blas import dgemm
import time

N = 1024
NITER = 30

A_array = random((N,N))
B_array = random((N,N))
C_array = zeros((N,N))

start = time.time()
for it in range(NITER):
	C_array = dot(A_array,B_array)
stop = time.time()

print"Time for 2d array dot product is:", (stop - start)/NITER

A_mat = matrix(A_array)
B_mat = matrix(B_array)
C_mat = matrix(C_array)

start = time.time()
for it in range(NITER):
	C_mat = A_mat*B_mat
stop = time.time()
print "Time for Matrix multiplication is:", (stop-start)/NITER

alpha = 1
start = time.time()
for it in range(NITER):
	C_array = dgemm(alpha, A_array, B_array)
stop = time.time()
print "Time for dgemm is:", (stop-start)/NITER

