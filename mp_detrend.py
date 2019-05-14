from __future__ import division
import untrendy
import george
#import cofiam
import astropy
import numpy as np
from cofiam import cofiam_iterative, max_order
from scipy.ndimage import median_filter



def cofiam_detrend(times, fluxes, errors, max_degree=30):
	best_model, best_degree, best_DW, max_degree = cofiam_iterative(times, fluxes, max_degree=max_degree)
	### detrend by dividing out the model
	flux_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return flux_detrend, errors_detrend 

def untrendy_detrend(times, fluxes, errors):
	f_detrend, sigma_detrend = untrendy.untrend(times, fluxes, errors)
	return f_detrend, sigma_detrend

def george_detrend(times, fluxes, errors):
	print('Nothing happening right now.')

def medfilt_detrend(fluxes, errors, size=5):
	detrend_errors = errors / fluxes
	detrend_fluxes = median_filter(fluxes, size=size, mode='nearest')
	return detrend_fluxes, detrend_errors

