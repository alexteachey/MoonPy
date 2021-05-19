from __future__ import division
import numpy as np
#from numba.core.decorators import jit
#from numba import float64, i8, jit
try:
	from numba.core.decorators import jit
except:
	from numba import jit 

"""
we need to solve the problem
AX = B

where A is a vector of coefficients for our linear problem
X is a matrix of terms, multiplying those values by the coefficients in A will give us
B the function values.



THIS CODE IS EXACTLY LIKE COFIAM, EXCEPT IT FITS POLYNOMIALS TO THE ENTIRE QUARTER INSTEAD OF SINUSOIDS.


NOTE THAT THE polynomial ALGORITHM FITS TERMS AS FOLLOWS

offset + (amp_s1 * (sin(2pi * time * 1) / (2 * baseline)) + amp_c1 * (cos(2*pi*time * 1) / 2*baseline) + ... up to the degree in question.

NOW FOR THE MATRIX REPRESENTATION, YOU NEED TO DO THIS FOR EVERY TIMESTEP! The matrix rows are for each time in your array!

"""

def max_order(times, duration, baseline=0, kmaximum=20):
	if baseline == 0:
		baseline = np.nanmax(times) - np.nanmin(times)
	kmax = int((2*baseline) / (12*duration))
	if kmax > kmaximum:
		kmax = kmaximum
	if kmax == 0:
		kmax = 1
	return kmax 



def DurbinWatson(residuals):
	residual_terms = []
	for nres, res in enumerate(residuals):
		try:
			residual_terms.append(residuals[nres+1] - residuals[nres])
		except:
			pass
	residual_terms = np.array(residual_terms)
	numerator = np.sum(residual_terms**2)
	denominator = np.sum(residuals**2)
	return numerator / denominator



def BIC(model, data, errors, nparams):
	chi2 = np.nansum(((model - data) / errors)**2)
	BICval = nparams*np.log(len(data)) + chi2
	return BICval


##########

# POLYNOMIAL
### AUTOCORRELATION
##### MINIMIZATION 

##########

### Special thanks to Michael Hippke for speeding this function up by orders of magnitude!
#@jit(fastmath=True, nopython=True, cache=True)
#@jit
#@jit((float64[:], i8))
@jit(debug=True, fastmath=True, nopython=True, cache=True)
def polyAM_matrix_gen(times, degree):
	baseline = np.nanmax(times) - np.nanmin(times)
	rows = len(times)
	#cols = 2 * (degree+1) #### this was the CoFiAM formulation 
	cols = (degree+1) #### if degree = 2, you need three columns -- for ax^2 + bx^1 + cx^0
	X_matrix = np.ones(shape=(rows,cols))
	for x in range(rows): ### for each row, that is, for each time!
		for y in range(1, int(cols/2)): #### for each column, up to half the columns -- why?
			#sinarg = (2*np.pi*times[x] * y) / baseline
			#X_matrix[x,y*2] = np.sin(sinarg)
			#X_matrix[x,y*2 + 1] = np.cos(sinarg)
			X_matrix[x,y*2] = times[x]**y #### TRY THIS!
		X_matrix[x,1] = times[x]
	return X_matrix 


def polyAM_matrix_coeffs(times, fluxes, degree):
	assert len(times) > 0
	Xmat = polyAM_matrix_gen(times, degree)
	beta_coefs = np.linalg.lstsq(Xmat, fluxes, rcond=None)[0]
	return Xmat, beta_coefs



### this function spits out the best fit line!
def polyAM_function(times, fluxes, degree):
	input_times = times.astype('f8')
	input_fluxes = fluxes.astype('f8')
	polyAM_matrix, polyAM_coefficients = polyAM_matrix_coeffs(input_times, input_fluxes, degree)
	output = np.matmul(polyAM_matrix, polyAM_coefficients)
	return output 


def polyAM_iterative(times, fluxes, max_degree=20, min_degree=1):
	### this function utilizes polyAM_function above, iterates it up to max_degree.
	### max degree may be calculated using max_order function

	vals_to_min = []
	degs_to_try = np.arange(min_degree,max_degree+1,1)
	DWstats = []

	for deg in degs_to_try:
		#print("k = ", deg)
		output_function = polyAM_function(times, fluxes, deg)

		residuals = fluxes - output_function

		DWstat = DurbinWatson(residuals)
		DWstats.append(DWstat)

		val_to_minimize = (DWstat - 2)**2
		vals_to_min.append(val_to_minimize)

	best_degree = degs_to_try[np.argmin(np.array(vals_to_min))]
	best_DW = DWstats[np.argmin(np.array(vals_to_min))]

	### re-generate the function with the best degree

	best_model = polyAM_function(times, fluxes, best_degree)

	return best_model, best_degree, best_DW, max_degree 


##########

# LOCAL
### POLYNOMIAL
##### FIT 

##########






### Special thanks to Michael Hippke for speeding this function up by orders of magnitude!
#@jit(fastmath=True, nopython=True, cache=True)
#@jit
#@jit((float64[:], i8))
@jit(debug=True, fastmath=True, nopython=True, cache=True)
def polyLOC_matrix_gen(times, degree):
	baseline = np.nanmax(times) - np.nanmin(times)
	rows = len(times)
	#cols = 2 * (degree+1) #### this was the CoFiAM formulation 
	cols = (degree+1) #### if degree = 2, you need three columns -- for ax^2 + bx^1 + cx^0
	X_matrix = np.ones(shape=(rows,cols))
	for x in range(rows): ### for each row, that is, for each time!
		for y in range(1, int(cols/2)): #### for each column, up to half the columns -- why?
			#sinarg = (2*np.pi*times[x] * y) / baseline
			#X_matrix[x,y*2] = np.sin(sinarg)
			#X_matrix[x,y*2 + 1] = np.cos(sinarg)
			X_matrix[x,y*2] = times[x]**y #### TRY THIS!
		X_matrix[x,1] = times[x]
	return X_matrix 


def polyLOC_matrix_coeffs(times, fluxes, degree):
	assert len(times) > 0
	Xmat = polyLOC_matrix_gen(times, degree)
	beta_coefs = np.linalg.lstsq(Xmat, fluxes, rcond=None)[0]
	return Xmat, beta_coefs



### this function spits out the best fit line!
def polyLOC_function(times, fluxes, degree):
	input_times = times.astype('f8')
	input_fluxes = fluxes.astype('f8')
	polyLOC_matrix, polyLOC_coefficients = polyLOC_matrix_coeffs(input_times, input_fluxes, degree)
	output = np.matmul(polyLOC_matrix, polyLOC_coefficients)
	return output 


def polyLOC_iterative(times, fluxes, errors, max_degree=20, min_degree=1):
	### this function utilizes polyAM_function above, iterates it up to max_degree.
	### max degree may be calculated using max_order function

	vals_to_min = []
	degs_to_try = np.arange(min_degree,max_degree+1,1)
	BICstats = []

	for deg in degs_to_try:
		output_function = polyAM_function(times, fluxes, deg) ### this is the model
		residuals = fluxes - output_function
		BICstat = BIC(output_function, fluxes, errors, deg+1)
		BICstats.append(BICstat)

	BICstats = np.array(BICstats)

	best_degree = degs_to_try[np.argmin(BICstats)]
	best_BIC = BICstats[np.argmin(np.array(BICstats))]

	### re-generate the function with the best degree

	best_model = polyLOC_function(times, fluxes, best_degree)

	return best_model, best_degree, best_BIC, max_degree 





