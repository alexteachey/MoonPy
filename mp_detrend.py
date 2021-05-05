from __future__ import division
import astropy
import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import medfilt
from scipy.interpolate import interp1d
from cofiam import cofiam_iterative
from poly_detrender import polyAM_iterative, polyLOC_iterative
import traceback
import matplotlib.pyplot as plt 


def cofiam_detrend(times, fluxes, errors, telescope='kepler', remove_outliers='y', outsig=3, window=19, mask_idxs=None, max_degree=30):
	print('calling mp_detrend.py/cofiam_detrend().')
	print("len(mask_idxs) [in-transit data] = ", len(mask_idxs))

	if type(mask_idxs) != type(None):
		if len(mask_idxs) > 0:
			unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		elif len(mask_idxs) == 0:
			##### ALL TIMES, FLUXES, AND ERRORS ARE UNMASKED.
			unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors

	if remove_outliers == 'y':
		outlier_idxs = []
		movmed = medfilt(unmasked_fluxes, kernel_size=window)
		for flidx, fl in enumerate(unmasked_fluxes):
			if np.abs(unmasked_fluxes[flidx] - movmed[flidx]) > outsig*unmasked_errors[flidx]:
				outlier_idxs.append(flidx)
		outlier_idxs = np.array(outlier_idxs)
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(unmasked_times, outlier_idxs), np.delete(unmasked_fluxes, outlier_idxs), np.delete(unmasked_errors, outlier_idxs)



	if telescope.lower() == 'tess':
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
		first_half_idxs = np.arange(0,largest_gap_idx+1,1)
		second_half_idxs = np.arange(largest_gap_idx+1,len(unmasked_times),1)

		try:
			if len(first_half_idxs) > 2:
				print('detrending first half of the quarter / section.')
				best_model1, best_degree1, best_DW1, max_degree1 = cofiam_iterative(np.array(unmasked_times[first_half_idxs], dtype=np.float64), np.array(unmasked_fluxes[first_half_idxs], dtype=np.float64), max_degree=int(max_degree))
			if len(second_half_idxs) > 2:
				print('detrending second half of the quarter / section.')
				best_model2, best_degree2, best_DW2, max_degree2 = cofiam_iterative(np.array(unmasked_times[second_half_idxs], dtype=np.float64), np.array(unmasked_fluxes[second_half_idxs], dtype=np.float64), max_degree=int(max_degree))
			
			if (len(first_half_idxs) > 2) and (len(second_half_idxs) > 2):
				best_model = np.concatenate((best_model1, best_model2))
				best_degree = np.nanmean((best_degree1, best_degree2))
				best_DW = np.nanmean((best_DW1, best_DW2))
				max_degree = np.nanmax((max_degree1, max_degree2))

			else:
				if len(first_half_idxs) < 2:
					best_model, best_degree, best_DW, max_degree = best_model2, best_degree2, best_DW2, max_degree2
					unmasked_times, unmasked_fluxes = unmasked_times[second_half_idxs], unmasked_fluxes[second_half_idxs]
				elif len(second_half_idxs) < 2:
					best_model, best_degree, best_DW, max_degree = best_model1, best_degree1, best_DW1, max_degree1		
					unmasked_times, unmasked_fluxes = unmasked_times[first_half_idxs], unmasked_fluxes[first_half_idxs]		


		except:
			traceback.print_exc()
			print('unable to call cofiam_iterative. Data points likely reduced to zero.')

	else: ### self.telescope != 'tess'
		try:
			best_model, best_degree, best_DW, max_degree = cofiam_iterative(np.array(unmasked_times, dtype=np.float64), np.array(unmasked_fluxes, dtype=np.float64), max_degree=int(max_degree))
			print(' ')
			print(' ')
		except:
			best_model, best_degree, best_DW, max_degree = np.array(unmasked_fluxes, dtype=np.float64), np.nan, np.nan, np.nan
			print('unable to call cofiam_iterative. Data points likely reduced to zero.')

	### at this point you have eliminated quite a few points, including the transit! So you need to interpolate to get the function values
	### at those locations in the time series, and to keep flux_detrend and errors_detrend the same length as the original time series.
	cofiam_interp = interp1d(unmasked_times, best_model, bounds_error=False, fill_value='extrapolate')
	try:
		best_model = cofiam_interp(times)
	except:
		best_model = cofiam_interp(np.array(times, dtype=np.float64))

	### detrend by dividing out the model
	flux_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return best_model, flux_detrend, errors_detrend #### best model is the COFIAM




def polyAM_detrend(times, fluxes, errors, telescope=None, remove_outliers='y', outsig=3, window=19, mask_idxs=None, max_degree=20):
	print('calling mp_detrend.py/polyAM_detrend().')
	if type(mask_idxs) != type(None):

		if len(mask_idxs) > 0:
			unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		
		else:
			unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors

	if remove_outliers == 'y':
		outlier_idxs = []
		movmed = medfilt(unmasked_fluxes, kernel_size=window)
		for flidx, fl in enumerate(unmasked_fluxes):
			if np.abs(unmasked_fluxes[flidx] - movmed[flidx]) > outsig*unmasked_errors[flidx]:
				outlier_idxs.append(flidx)
		outlier_idxs = np.array(outlier_idxs)
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(unmasked_times, outlier_idxs), np.delete(unmasked_fluxes, outlier_idxs), np.delete(unmasked_errors, outlier_idxs)


	if telescope.lower() == 'tess':
		### you need to detrend the two halves separately!
		deltats = []
		for nut, ut in enumerate(unmasked_times):
			if nut == unmasked_times.size-1:
				pass
			else:
				deltats.append(unmasked_times[nut+1] - unmasked_times[nut])

		deltats = np.array(deltats)
		deltat = np.nanmedian(deltats)


		### identify the largest gap! -- THIS SHOULD BE THE DATA GAP IN THE TESS SECTOR.
		largest_gap_idx = np.nanargmax(deltats) ### this will be the index of the last(!) data point before the gap
		first_half_idxs = np.arange(0,largest_gap_idx+1,1)
		second_half_idxs = np.arange(largest_gap_idx+1,len(unmasked_times),1)

		try:
			if len(first_half_idxs) > 2:
				print('detrending first section.')
				best_model1, best_degree1, best_DW1, max_degree1 = polyAM_iterative(np.array(unmasked_times[first_half_idxs], dtype=np.float64), np.array(unmasked_fluxes[first_half_idxs], dtype=np.float64), max_degree=int(max_degree))
			
			if len(second_half_idxs) > 2:
				print('detrending second section.')
				best_model2, best_degree2, best_DW2, max_degree2 = polyAM_iterative(np.array(unmasked_times[second_half_idxs], dtype=np.float64), np.array(unmasked_fluxes[second_half_idxs], dtype=np.float64), max_degree=int(max_degree))
			
			if (len(first_half_idxs) > 2) and (len(second_half_idxs) > 2):
				best_model = np.concatenate((best_model1, best_model2))
				best_degree = np.nanmean((best_degree1, best_degree2))
				best_DW = np.nanmean((best_DW1, best_DW2))
				max_degree = np.nanmax((max_degree1, max_degree2))

			else:
				if len(first_half_idxs) < 2:
					best_model, best_degree, best_DW, max_degree = best_model2, best_degree2, best_DW2, max_degree2
					unmasked_times, unmasked_fluxes = unmasked_times[second_half_idxs], unmasked_fluxes[second_half_idxs]
				elif len(second_half_idxs) < 2:
					best_model, best_degree, best_DW, max_degree = best_model1, best_degree1, best_DW1, max_degree1		
					unmasked_times, unmasked_fluxes = unmasked_times[first_half_idxs], unmasked_fluxes[first_half_idxs]						



		except:
			traceback.print_exc()
			print('unable to call polyAM_iterative. Data points likely reduced to zero.')


	else:
		try:
			best_model, best_degree, best_DW, max_degree = polyAM_iterative(np.array(unmasked_times, dtype=np.float64), np.array(unmasked_fluxes, dtype=np.float64), max_degree=int(max_degree))
			print(' ')
			print(' ')
		except:
			traceback.print_exc()
			print('unable to call polyAM_iterative. Data points likely reduced to zero.')



	### at this point you have eliminated quite a few points, including the transit! So you need to interpolate to get the function values
	### at those locations in the time series, and to keep flux_detrend and errors_detrend the same length as the original time series.
	polyAM_interp = interp1d(unmasked_times, best_model, bounds_error=False, fill_value='extrapolate')
	try:
		best_model = polyAM_interp(times)
	except:
		best_model = polyAM_interp(np.array(times, dtype=np.float64))

	### detrend by dividing out the model
	flux_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return best_model, flux_detrend, errors_detrend 





def polyLOC_detrend(times, fluxes, errors, telescope=None, remove_outliers='y', outsig=3, window=19, local_window_duration_multiple=10, Tmid=None, Tdur_days=None, mask_idxs=None, max_degree=20):
	print('calling mp_detrend.py/polyLOC_detrend().')
	#### THIS DETRENDING HAS TO BE USED ON A TRANSIT BY TRANSIT BASIS, *NOT* A QUARTER BY QUARTER BASIS.

	#### use the local_window_duration_multiple to use only times from transit midtime
	if (Tdur_days != None) and (Tmid != None):
		#### you need these to calculate
		local_window_days = Tdur_days * local_window_duration_multiple
		local_window_nobs = local_window_days * 48 #### roughly 48 observations per day (once every half hour)
		Tmid_idx = np.nanargmin(np.abs(times - Tmid))
		window_idxs = np.arange(Tmid_idx - local_window_nobs, Tmid_idx + local_window_nobs, 1)
		times, fluxes, errors = times[window_idxs], fluxes[window_idxs], errors[window_idxs]

	if type(mask_idxs) != type(None):

		if len(mask_idxs) > 0:
			unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		else:
			unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors

	if remove_outliers == 'y':
		outlier_idxs = []
		movmed = medfilt(unmasked_fluxes, kernel_size=window)
		for flidx, fl in enumerate(unmasked_fluxes):
			if np.abs(unmasked_fluxes[flidx] - movmed[flidx]) > outsig*unmasked_errors[flidx]:
				outlier_idxs.append(flidx)
		outlier_idxs = np.array(outlier_idxs)
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(unmasked_times, outlier_idxs), np.delete(unmasked_fluxes, outlier_idxs), np.delete(unmasked_errors, outlier_idxs)

	if telescope.lower() == 'tess':
		### you need to detrend the two halves separately!
		deltats = []
		for nut, ut in enumerate(unmasked_times):
			if nut == unmasked_times.size-1:
				pass
			else:
				deltats.append(unmasked_times[nut+1] - unmasked_times[nut])
		deltats = np.array(deltats)
		deltat = np.nanmedian(deltats)


		### identify the largest gap! -- THIS SHOULD BE THE DATA GAP IN THE TESS SECTOR.
		largest_gap_idx = np.nanargmax(deltats) ### this will be the index of the last(!) data point before the gap
		first_half_idxs = np.arange(0,largest_gap_idx+1,1)
		second_half_idxs = np.arange(largest_gap_idx+1,len(unmasked_times),1)

		try:
			best_model1, best_degree1, best_BIC1, max_degree1 = polyLOC_iterative(np.array(unmasked_times[first_half_idxs], dtype=np.float64), np.array(unmasked_fluxes[first_half_idxs], dtype=np.float64), max_degree=int(max_degree))
			best_model2, best_degree2, best_BIC2, max_degree2 = polyLOC_iterative(np.array(unmasked_times[second_half_idxs], dtype=np.float64), np.array(unmasked_fluxes[second_half_idxs], dtype=np.float64), max_degree=int(max_degree))
			best_model = np.concatenate((best_model1, best_model2))
			best_degree = np.nanmean((best_degree1, best_degree2))
			best_BIC = np.nanmean((best_BIC1, best_BIC2))
			max_degree = np.nanmax((max_degree1, max_degree2))
		except:
			traceback.print_exc()
			print('unable to call polyAM_iterative. Data points likely reduced to zero.')

	else:
		try:
			best_model, best_degree, best_BIC, max_degree = polyLOC_iterative(np.array(unmasked_times, dtype=np.float64), np.array(unmasked_fluxes, dtype=np.float64), np.array(unmasked_errors, dtype=np.float64), max_degree=int(max_degree))
		except:
			traceback.print_exc()
			print('unable to call polyAM_iterative. Data points likely reduced to zero.')

	### at this point you have eliminated quite a few points, including the transit! So you need to interpolate to get the function values
	### at those locations in the time series, and to keep flux_detrend and errors_detrend the same length as the original time series.
	polyLOC_interp = interp1d(unmasked_times, best_model, bounds_error=False, fill_value='extrapolate')
	try:
		best_model = polyLOC_interp(times)
	except:
		best_model = polyLOC_interp(np.array(times, dtype=np.float64))

	### detrend by dividing out the model
	flux_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return best_model, times, flux_detrend, errors_detrend 





def untrendy_detrend(times, fluxes, errors, telescope=None, mask_idxs=None):
	print('calling mp_detrend.py/untrendy_detrend().')
	import untrendy

	print('BEWARE: Untrendy is failing because of a strange bug within scipy.')
	if type(mask_idxs) != type(None):
		print(' ')
		print(type(mask_idxs))
		print(' ')

		if (type(mask_idxs) != type(None)) and (len(mask_idxs) > 0):
			unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		else:
			unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors
		### untrendy throws an error if unmasked_times arent strictly increasing
		time_diffs = np.diff(unmasked_times)
		if np.any(time_diffs <= 0):
			print("TIMES ARE NOT STRICTLY INCREASING!")
		f_detrend = untrendy.untrend(unmasked_times, unmasked_fluxes, yerr=unmasked_errors)[0]
		untrendy_interp = interp1d(unmasked_times, f_detrend, bounds_error=False, fill_value='extrapolate')
		f_detrend = untrendy_interp(times)
		sigma_detrend = untrendy.untrend(times, fluxes, errors)[1]
		best_model = untrendy.fit_trend(times, fluxes, errors)
	else:
		f_detrend, sigma_detrend = untrendy.untrend(times, fluxes, errors)
		best_model = untrendy.fit_trend(times, fluxes, errors)
	return best_model, f_detrend, sigma_detrend





def george_detrend(times, fluxes, errors, GP_kernel='ExpSquaredKernel', metric=1.0, telescope=None, mask_idxs=None):
	print('calling mp_detrend.py/george_detrend().')
	import george

	if GP_kernel != 'ExpSquaredKernel':
		try:
			kernel_choice = vars(george.kernels)[GP_kernel] ### accesses the kernel through a dictionary, with kernel_name being the key.
			print('using ', GP_kernel)
		
		except:
			print('GP_kernel input in self.detrend() or george_detrend() is missing or invalid. unable to load your kernel choice. Loading ExpSquaredKernel.')
			from george.kernels import ExpSquaredKernel as kernel_choice
			print("george GP code is using the Exponential Squared Kernel with metric="+str(metric)+'.')
		
	elif GP_kernel == 'ExpSquaredKernel':
		from george.kernels import ExpSquaredKernel as kernel_choice
		print("george GP code is using the Exponential Squared Kernel with metric="+str(metric))


	elif GP_kernel == None:
		from george.kernels import ExpSquaredKernel as kernel_choice
		print("george GP code is using the Exponential Squared Kernel with metric="+str(metric)+'.')

	print('mask_idxs = ', mask_idxs)
	if (type(mask_idxs) != type(None)) and (len(mask_idxs) > 0):
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
	else:
		unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors 

	kernel_arg = np.var(unmasked_fluxes) * kernel_choice(metric=metric)
	print('generating the gp...')
	gp = george.GP(kernel_arg)
	print('computing the gp...')
	gp.compute(unmasked_times, unmasked_errors) ### pre-compute the factorization of the matrix

	### compute the log likelihood
	print('computing gp.lnlikelihood...')
	print(gp.lnlikelihood(unmasked_fluxes))

	### now interpolate
	print('predicting...')
	gp_mu, gp_cov = gp.predict(unmasked_fluxes, times)  ### FIRST ARGUMENT ARE THE ORIGINAL y-values, *NOT* ALL y-values! (you leave out the transit times)
	gp_std = np.sqrt(np.diag(gp_cov))
	flux_detrend = fluxes / gp_mu 
	errors_detrend = errors / fluxes 
	best_model = gp_mu
	print(' ')
	print(' ')

	return best_model, flux_detrend, errors_detrend





def medfilt_detrend(times, fluxes, errors, kernel_hours, telescope=None, mask_idxs=None):
	print('calling np_detrend.py/medfilt_detrend().')

	print('kernel_hours = ', kernel_hours)

	kernel_size = int(2*kernel_hours) #### 2 data points per hour for Kepler.
	if kernel_size % 2 == 0:
		kernel_size = kernel_size + 1
		assert kernel_size % 2 == 1

	if (type(mask_idxs) != type(None)) and (len(mask_idxs) > 0):
		print('performing median filter with masked points.')
		### that is, if there are masks for the transits (there should be!)
		unmasked_times, unmasked_fluxes, unmasked_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
	else:
		unmasked_times, unmasked_fluxes, unmasked_errors = times, fluxes, errors

		print('performing median filter without masked points.')

	try:
		flux_trend = median_filter(unmasked_fluxes, size=kernel_size, mode='nearest')
		print('utilizing scipy.ndimage.median_filter().')
	
	except:
		flux_trend = medfilt(unmasked_fluxes, kernel_size=kernel_size)
		print('median_filter() failed, utilizing scipy.signal.medfilt().')

	medfilt_interp = interp1d(unmasked_times, flux_trend, bounds_error=False, fill_value='extrapolate')
	flux_trend = medfilt_interp(times)	
	best_model = flux_trend

	detrend_errors = errors / fluxes 
	detrend_fluxes = fluxes / flux_trend

	return best_model, detrend_fluxes, detrend_errors





def methmarg_detrend(times, fluxes, errors, kernel_hours, GP_kernel='ExpSquaredKernel', metric=1.0, local_window_duration_multiple=10, Tdur_days=None, Tmid=None, telescope=None, mask_idxs=None, max_degree=30):
	print('calling np_detrend.py/methmarg_detrend().')
	#### THIS FUNCTION WILL RUN ALL THE OTHER DETRENDERS, and take the median value at every time step.
	try:
		cofiam_fluxes, cofiam_errors = cofiam_detrend(times=times, fluxes=fluxes, errors=errors, telescope=telescope, remove_outliers='y', outsig=3, window=19, mask_idxs=mask_idxs, max_degree=30)
		include_cofiam = 'y'
	except:
		traceback.print_exc()
		include_cofiam = 'n'

	try:
		polyAM_fluxes, polyAM_errors = polyAM_detrend(times=times, fluxes=fluxes, errors=errors, telescope=telescope, remove_outliers='y', outsig=3, window=19, mask_idxs=mask_idxs, max_degree=20)	
		include_polyAM = 'y'
	except:
		traceback.print_exc()
		include_polyAM = 'n'

	try:
		polyLOC_times, polyLOC_fluxes, polyLOC_errors = polyLOC_detrend(times=times, fluxes=fluxes, errors=errors, telescope=telescope, remove_outliers='y', outsig=3, window=19, local_window_duration_multiple=local_window_duration_multiple, Tmid=Tmid, Tdur_days=Tdur_days, mask_idxs=mask_idxs, max_degree=20)
		include_polyLOC = 'y'

		#### need to embed this in a larger array, so that it's the same size.
		full_polyLOC_fluxes, full_polyLOC_errors = np.full(len(times), np.nan), np.full(len(times), np.nan)
		match_idxs = []
		for i in np.arange(0,len(times),1):
			if times[i] in polyLOC_times:
				match_idx = np.where(times[i] == polyLOC_times)[0]
				full_polyLOC_fluxes[i] = polyLOC_fluxes[match_idx]
				full_polyLOC_errors[i] = polyLOC_errors[match_idx]

		polyLOC_fluxes, polyLOC_errors = full_polyLOC_fluxes, full_polyLOC_errors #### should be the same size as the other arrays!
		
	except:
		traceback.print_exc()
		include_polyLOC = 'n'

	try:
		medfilt_fluxes, medfilt_errors = medfilt_detrend(times=times, fluxes=fluxes, errors=errors, kernel_hours=kernel_hours, telescope=telescope, mask_idxs=mask_idxs)	
		include_medfilt = 'y'
	except:
		traceback.print_exc()
		include_medfilt = 'n'

	methmarg_fluxes = []
	mathmarg_errors = []


	#### IDEALLY THIS DETRENDING IS ALREADY HAPPENING ON A QUARTER BY QUARTER BASIS -- YOU ARE NOT 
	if include_cofiam == 'y':
		quarter_flux_stack = cofiam_fluxes
		quarter_error_stack = cofiam_errors
	else:
		pass

	if (include_polyAM == 'y') and (include_cofiam == 'y'):
		quarter_flux_stack = np.vstack((quarter_flux_stack, polyAM_fluxes))
		quarter_error_stack = np.vstack((quarter_error_stack, polyAM_errors))
	elif (include_polyAM == 'y') and (include_cofiam == 'n'):
		quarter_flux_stack = polyAM_fluxes 
		quarter_error_stack = polyAM_errors
	elif include_polyAM == 'n':
		pass

	if (include_polyLOC == 'y') and ((include_cofiam == 'y') or (include_polyAM == 'y')):
		quarter_flux_stack = np.vstack((quarter_flux_stack, polyLOC_fluxes))
		quarter_error_stack = np.vstack((quarter_error_stack, polyLOC_errors))
	elif (include_polyLOC == 'y') and include_cofiam == 'n' and include_polyAM == 'n':
		quarter_flux_stack = polyLOC_fluxes 
		quarter_error_stack = polyLC_errors 
	elif include_polyLOC == 'n':
		pass

	if (include_medfilt == 'y') and ((include_cofiam == 'y') or (include_polyAM == 'y') or (include_polyLOC == 'y')):
		quarter_flux_stack = np.vstack((quarter_flux_stack, medfilt_fluxes))
		quarter_error_stack = np.vstack((quarter_error_stack, medfilt_errors))
	elif (include_medfilt == 'y') and (include_cofiam == 'n') and (include_polyAM == 'n') and (include_polyLOC == 'n'):
		quarter_flux_stack = medfilt_fluxes 
		quarter_error_stack = medfilt_errors
	else:
		pass

	print("number of methods included in MethMarg = ", quarter_flux_stack.shape[0])

	quarter_methmarg_median_fluxes = np.nanmedian(quarter_flux_stack, axis=0)
	quarter_methmarg_median_errors = np.zeros(shape=quarter_methmarg_median_fluxes.shape)
	for i in np.arange(0,quarter_flux_stack.shape[1],1):
		iMAD = np.nanmedian(np.abs(quarter_error_stack.T[i] - np.nanmedian(quarter_error_stack.T[i])))
		isig = 1.4826 * iMAD
		quarter_methmarg_median_errors[i] = isig 


	print('quarter_methmarg_median_fluxes.shape = ', quarter_methmarg_median_fluxes.shape)
	print('quarter_methmarg_median_errors.shape = ', quarter_methmarg_median_errors.shape)

	#### have to duplicate the fluxes here to take the place of "best_model" output from others
	return quarter_methmarg_median_fluxes, quarter_methmarg_median_fluxes, quarter_methmarg_median_errors


