from __future__ import division
import astropy
import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import medfilt
from scipy.interpolate import interp1d
from cofiam import cofiam_iterative, max_order



def cofiam_detrend(times, fluxes, errors, remove_outliers='y', outsig=3, window=19, mask_idxs=None, max_degree=30):
	if type(mask_idxs) != type(None):
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)

		if remove_outliers == 'y':
			outlier_idxs = []
			movmed = medfilt(unmasked_fluxes, kernel_size=window)
			for flidx, fl in enumerate(unmasked_fluxes):
				if np.abs(unmasked_fluxes[flidx] - movmed[flidx]) > outsig*unmasked_errors[flidx]:
					outlier_idxs.append(flidx)
			outlier_idxs = np.array(outlier_idxs)
			unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(unmasked_times, outlier_idxs), np.delete(unmasked_fluxes, outlier_idxs), np.delete(unmasked_errors, outlier_idxs)

		try:
			best_model, best_degree, best_DW, max_degree = cofiam_iterative(unmasked_times, unmasked_fluxes, max_degree=max_degree)
		except:
			print('unable to call cofiam_iterative. Data points likely reduced to zero.')
		### at this point you have eliminated quite a few points, including the transit! So you need to interpolate to get the function values
		### at those locations in the time series, and to keep flux_detrend and errors_detrend the same length as the original time series.
		cofiam_interp = interp1d(unmasked_times, best_model, bounds_error=False, fill_value='extrapolate')
		best_model = cofiam_interp(times)
	else:
		best_model, best_degree, best_DW, max_degree = cofiam_iterative(times, fluxes, max_degree=max_degree)
	### detrend by dividing out the model
	flux_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return flux_detrend, errors_detrend 

def untrendy_detrend(times, fluxes, errors, mask_idxs=None):
	import untrendy

	if type(mask_idxs) != type(None):
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		f_detrend = untrendy.detrend(unmasked_times, unmasked_fluxes, unmasked_errors)[0]
		untrendy_interp = interp1d(unmasked_times, f_detrend, bounds_error=False, fill_value='extrapolate')
		f_detrend = untrendy_interp(times)
		sigma_detrend = untrendy.untrend(times, fluxes, errors)[1]
	else:
		f_detrend, sigma_detrend = untrendy.untrend(times, fluxes, errors)
	return f_detrend, sigma_detrend

def george_detrend(times, fluxes, errors, mask_idxs=None):
	import george

	print('Nothing happening right now.')

def medfilt_detrend(times, fluxes, errors, size, mask_idxs=None):
	if type(mask_idxs) != type(None):
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		flux_trend = median_filter(fluxes, size=size, mode='nearest')
		medfilt_interp = interp1d(unmasked_times, unmasked_fluxes, bounds_error=False, fill_value='extrapolate')
		flux_trend = medfilt_interp(times)

	if size == None:
		size = 9 ### standardize it
	detrend_errors = errors / fluxes
	#flux_trend = median_filter(fluxes, size=size, mode='nearest')
	flux_trend = medfilt(fluxes, kernel=size)
	detrend_fluxes = fluxes / flux_trend

	return detrend_fluxes, detrend_errors

