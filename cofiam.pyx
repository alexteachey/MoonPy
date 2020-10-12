from __future__ import division
import numpy as np
from numba import jit

"""
we need to solve the problem
AX = B

where A is a vector of coefficients for our linear problem
X is a matrix of terms, multiplying those values by the coefficients in A will give us
B the function values.


NOTE THAT THE COFIAM ALGORITHM FITS TERMS AS FOLLOWS

offset + (amp_s1 * (sin(2pi * time * 1) / (2 * baseline)) + amp_c1 * (cos(2*pi*time * 1) / 2*baseline) + ... up to the degree in question.

NOW FOR THE MATRIX REPRESENTATION, YOU NEED TO DO THIS FOR EVERY TIMESTEP! The matrix rows are for each time in your array!

"""

def max_order(times, duration, baseline=0, kmaximum=30):
	#cdef float[:] ctimes = times
	ctimes = times 
	#cdef double[:] cduration = duration
	cduration = duration
	cdef int ckmaximum = kmaximum

	if baseline == 0:
		baseline = np.nanmax(ctimes) - np.nanmin(ctimes)
	kmax = int((2*baseline) / (12*duration))
	if kmax > ckmaximum:
		kmax = ckmaximum
	if kmax == 0:
		kmax = 1
	return kmax 



def DurbinWatson(residuals):
	cdef double[:] cresiduals = residuals

	residual_terms = []
	for nres, res in enumerate(cresiduals):
		try:
			residual_terms.append(cresiduals[nres+1] - cresiduals[nres])
		except:
			pass
	residual_terms = np.array(residual_terms)
	numerator = np.sum(residual_terms**2)
	denominator = np.sum(residuals**2)
	return numerator / denominator



### Special thanks to Michael Hippke for speeding this function up by orders of magnitude!
#### UPDATE -- October 8th -- @JIT is causing segmentation fault during detrending, no idea why. Will need to solve this.
#@jit(fastmath=True, nopython=True, cache=True)
#@jit
def cofiam_matrix_gen(times, degree):

	cdef double[:] ctimes = times
	cdef int cdegree = degree

	baseline = np.nanmax(ctimes) - np.nanmin(ctimes)
	rows = len(ctimes)
	cols = 2 * (cdegree+1)
	X_matrix = np.ones(shape=(rows,cols))
	for x in range(rows):
		for y in range(1, int(cols/2)):
			sinarg = (2*np.pi*ctimes[x] * y) / baseline
			X_matrix[x,y*2] = np.sin(sinarg)
			X_matrix[x,y*2 + 1] = np.cos(sinarg)
		X_matrix[x,1] = times[x]
	return X_matrix 


def cofiam_matrix_coeffs(times, fluxes, degree):
	cdef double[:] ctimes = times
	cdef double[:] cfluxes = fluxes
	cdef int cdegree = degree

	assert len(ctimes) > 0
	Xmat = cofiam_matrix_gen(ctimes, cdegree)
	beta_coefs = np.linalg.lstsq(Xmat, cfluxes, rcond=None)[0]
	return Xmat, beta_coefs



### this function spits out the best fit line!
def cofiam_function(times, fluxes, degree):
	cdef double[:] input_times = times
	cdef double[:] input_fluxes = fluxes
	cdef int cdegree = degree

	#input_times = times.astype('f8')
	#input_fluxes = fluxes.astype('f8')
	cofiam_matrix, cofiam_coefficients = cofiam_matrix_coeffs(input_times, input_fluxes, degree)
	output = np.matmul(cofiam_matrix, cofiam_coefficients)
	return output 


def cofiam_iterative(times, fluxes, max_degree=30, min_degree=1):
	### this function utilizes cofiam_function above, iterates it up to max_degree.
	### max degree may be calculated using max_order function

	cdef double[:] ctimes = times
	cdef double[:] cfluxes = fluxes
	cdef int cmax_degree = max_degree
	cdef int cmin_degree = min_degree

	vals_to_min = []
	degs_to_try = np.arange(cmin_degree,cmax_degree+1,1)
	DWstats = []

	for deg in degs_to_try:
		print("k = ", deg)
		output_function = cofiam_function(ctimes, cfluxes, deg)

		residuals = cfluxes - output_function

		DWstat = DurbinWatson(residuals)
		DWstats.append(DWstat)

		val_to_minimize = (DWstat - 2)**2
		vals_to_min.append(val_to_minimize)

	best_degree = degs_to_try[np.argmin(np.array(vals_to_min))]
	best_DW = DWstats[np.argmin(np.array(vals_to_min))]

	cdef int cbest_degree = best_degree

	### re-generate the function with the best degree

	best_model = cofiam_function(ctimes, cfluxes, cbest_degree)

	return best_model, best_degree, best_DW, max_degree 






