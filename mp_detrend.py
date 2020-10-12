from __future__ import division
import astropy
import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import medfilt
from scipy.interpolate import interp1d
from cofiam import cofiam_iterative, max_order
import traceback
import matplotlib.pyplot as plt 
import pyximport 




def cofiam_detrend(times, fluxes, errors, telescope=None, remove_outliers='y', outsig=3, window=19, mask_idxs=None, max_degree=30):
	print("len(mask_idxs) = ", len(mask_idxs))
	print('len(times) = ', len(times))

	if type(mask_idxs) != type(None):

		if len(mask_idxs) > 0:
			mask_idxs = np.array(mask_idxs, dtype=np.int8)
			print('type(mask_idxs) = ', type(mask_idxs))
			unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		else:
			print('mask_idxs = ', mask_idxs)
			unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors

	if remove_outliers == 'y':
		outlier_idxs = []
		movmed = medfilt(unmasked_fluxes, kernel_size=window)
		for flidx, fl in enumerate(unmasked_fluxes):
			if np.abs(unmasked_fluxes[flidx] - movmed[flidx]) > outsig*unmasked_errors[flidx]:
				outlier_idxs.append(flidx)
		outlier_idxs = np.array(outlier_idxs)
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(unmasked_times, outlier_idxs), np.delete(unmasked_fluxes, outlier_idxs), np.delete(unmasked_errors, outlier_idxs)

	if (telescope == 'tess') or (telescope == 'TESS') or (telescope == 'Tess'):
		### you need to detrend the two halves separately!
		deltats = []
		for nut, ut in enumerate(unmasked_times):
			if nut == unmasked_times.size-1:
				pass
			else:
				deltats.append(unmasked_times[nut+1] - unmasked_times[nut])

		deltats = np.array(deltats)
		deltat = np.nanmedian(deltats)

		### identify the largest gap!
		largest_gap_idx = np.nanargmax(deltats) ### this will be the index of the last(!) data point before the gap
		print('largest_gap_idx = ', largest_gap_idx)
		first_half_idxs = np.arange(0,largest_gap_idx+1,1)
		second_half_idxs = np.arange(largest_gap_idx+1,len(unmasked_times),1)

		try:
			best_model1, best_degree1, best_DW1, max_degree1 = cofiam_iterative(unmasked_times[first_half_idxs], unmasked_fluxes[first_half_idxs], max_degree=max_degree)		
			best_model2, best_degree2, best_DW2, max_degree2 = cofiam_iterative(unmasked_times[second_half_idxs], unmasked_fluxes[second_half_idxs], max_degree=max_degree)
			best_model = np.concatenate((best_model1, best_model2))
			best_degree = np.nanmean((best_degree1, best_degree2))
			best_DW = np.nanmean((best_DW1, best_DW2))
			max_degree = np.nanmax((max_degree1, max_degree2))
		except:
			traceback.print_exc()
			print('unable to call cofiam_iterative. Data points likely reduced to zero.')

	else:
		try:
			best_model, best_degree, best_DW, max_degree = cofiam_iterative(unmasked_times, unmasked_fluxes, max_degree=max_degree)
			print(' ')
			print(' ')
		except:
			traceback.print_exc()
			print('unable to call cofiam_iterative. Data points likely reduced to zero.')

	### at this point you have eliminated quite a few points, including the transit! So you need to interpolate to get the function values
	### at those locations in the time series, and to keep flux_detrend and errors_detrend the same length as the original time series.
	cofiam_interp = interp1d(unmasked_times, best_model, bounds_error=False, fill_value='extrapolate')
	best_model = cofiam_interp(times)

	### detrend by dividing out the model
	flux_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return flux_detrend, errors_detrend 




def untrendy_detrend(times, fluxes, errors, telescope=None, mask_idxs=None):
	import untrendy

	print('BEWARE: Untrendy is failing because of a strange bug within scipy.')
	if type(mask_idxs) != type(None):
		mask_idxs = np.array(mask_idxs, dtype=np.int8)
		print(' ')
		print(type(mask_idxs))
		print(' ')

		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		### untrendy throws an error if unmasked_times arent strictly increasing
		time_diffs = np.diff(unmasked_times)
		if np.any(time_diffs <= 0):
			print("TIMES ARE NOT STRICTLY INCREASING!")
		f_detrend = untrendy.untrend(unmasked_times, unmasked_fluxes, yerr=unmasked_errors)[0]
		untrendy_interp = interp1d(unmasked_times, f_detrend, bounds_error=False, fill_value='extrapolate')
		f_detrend = untrendy_interp(times)
		sigma_detrend = untrendy.untrend(times, fluxes, errors)[1]
	else:
		f_detrend, sigma_detrend = untrendy.untrend(times, fluxes, errors)
	return f_detrend, sigma_detrend

def george_detrend(times, fluxes, errors, telescope=None, mask_idxs=None):
	import george
	from george.kernels import ExpSquaredKernel

	unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)

	kernel = ExpSquaredKernel(1.0)
	gp = george.GP(kernel)
	gp.compute(unmasked_times, unmasked_errors) ### pre-compute the factorization of the matrix

	### compute the log likelihood
	print(gp.lnlikelihood(unmasked_fluxes))

	### now interpolate
	gp_mu, gp_cov = gp.predict(fluxes, times)
	gp_std = np.sqrt(np.diag(gp_cov))
	flux_detrend = fluxes / gp_mu 
	errors_detrend = errors / fluxes 

	return flux_detrend, errors_detrend



def medfilt_detrend(times, fluxes, errors, size, telescope=None, mask_idxs=None):
	if size == None:
		size = 9 

	print('kernel size = ', size)
	if type(mask_idxs) != type(None):
		print('performing median filter with masked points.')
		### that is, if there are masks for the transits (there should be!)
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		print('len(unmasked_times), len(unmasked_fluxes) = ', len(unmasked_times), len(unmasked_fluxes))
		try:
			flux_trend = median_filter(unmasked_fluxes, size=size, mode='nearest')
			print('utilizing scipy.ndimage.median_filter().')
		except:
			flux_trend = medfilt(unmasked_fluxes, kernel_size=size)
		print('median_filter() failed, utilizing scipy.signal.medfilt().')
		medfilt_interp = interp1d(unmasked_times, flux_trend, bounds_error=False, fill_value='extrapolate')
		flux_trend = medfilt_interp(times)

	else:
		print('performing median filter without masked points.')
		try:
			flux_trend = median_filter(fluxes, size=size, mode='nearest')
		except:
			flux_trend = medfilt(fluxes, kernel_size=size)

	detrend_errors = errors / fluxes 
	detrend_fluxes = fluxes / flux_trend

	return detrend_fluxes, detrend_errors

