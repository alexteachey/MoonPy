from __future__ import division
import numpy as np

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
	if baseline == 0:
		baseline = np.nanmax(times) - np.nanmin(times)
	kmax = int((2*baseline) / (12*duration))
	if kmax > kmaximum:
		kmax = kmaximum
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




def cofiam_matrix_gen(times, degree):
	baseline = np.nanmax(times) - np.nanmin(times)

	for ntval, tval in enumerate(times):
		matrix_row = [1, tval]
		for deg in np.arange(1,degree+1,1):
			sinterm = np.sin((2*np.pi*tval*deg)/(2*baseline))
			costerm = np.cos((2*np.pi*tval*deg)/(2*baseline))
			matrix_row.append(sinterm)
			matrix_row.append(costerm)
		matrix_row = np.array(matrix_row)

		if ntval == 0:
			X_matrix = matrix_row
		else:
			X_matrix = np.vstack((X_matrix, matrix_row))

	return X_matrix



def cofiam_matrix_coeffs(times, fluxes, degree):
	Xmat = cofiam_matrix_gen(times, degree)
	beta_coefs = np.linalg.lstsq(Xmat, fluxes)[0]
	return beta_coefs



### this function spits out the best fit line!
def cofiam_function(times, fluxes, degree):
	cofiam_matrix = cofiam_matrix_gen(times, degree)
	cofiam_coefficients = cofiam_matrix_coeffs(times, fluxes, degree)
	output = np.matmul(cofiam_matrix, cofiam_coefficients)
	return output 


def cofiam_iterative(times, fluxes, max_degree=30, min_degree=1):
	### this function utilizes cofiam_function above, iterates it up to max_degree.
	### max degree may be calculated using max_order function

	vals_to_min = []
	degs_to_try = np.arange(min_degree,max_degree+1,1)
	DWstats = []

	for deg in degs_to_try:
		print("k = ", deg)
		output_function = cofiam_function(times, fluxes, deg)

		residuals = fluxes - output_function

		DWstat = DurbinWatson(residuals)
		DWstats.append(DWstat)

		val_to_minimize = (DWstat - 2)**2
		vals_to_min.append(val_to_minimize)

	best_degree = degs_to_try[np.argmin(np.array(vals_to_min))]
	best_DW = DWstats[np.argmin(np.array(vals_to_min))]

	### re-generate the function with the best degree

	best_model = cofiam_function(times, fluxes, best_degree)

	return best_model, best_degree, best_DW, max_degree 






