""" Remember the Python zen: There should be only one obvious way to do 
things. We are going to compare the low-level LAPACK functions 
degtrf and dgetrs to the built in scipy LU routines. """

from numpy import array
from scipy.linalg import lu_solve, lu_factor
from scipy.linalg.lapack import dgetrf, dgetrs
from numpy.random import random
import time

# Make arrays
NUM_ITER = 10
N = 1024
A = random((N,N))
b = random((N,1)) 

# Solve using scipy.linalg
start = time.time()
for it in range(NUM_ITER):
	(LU_and_piv) = lu_factor(A)
	(x) = lu_solve(LU_and_piv, b)
stop = time.time()
print "Time for scipy routine is", (stop - start)/NUM_ITER

# Solve using scipy.linalg.lapack
start = time.time()
for it in range(NUM_ITER):
	(LU, piv, info) = dgetrf(A)
	(x) = dgetrs(LU, piv, b)
stop = time.time()
print "Time for LAPACK routine is", (stop - start)/NUM_ITER
