from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import pymultinest
import emcee
import untrendy
import sklearn
import kplr
import everest
import eleanor
import corner
#import forecaster
#import pyluna
#import tensorflow
#import cofiam
import astropy
import warnings

### The packages below interface with standard packages, but within the context
### of what you want here... maybe change this usage below!
from mp_lcfind import kplr_target_download, kplr_coord_download, eleanor_target_download, eleanor_coord_download
from mp_detrend import untrendy_detrend, cofiam_detrend, george_detrend
from mp_fit import mp_multinest, mp_emcee
#from pyluna import run_LUNA
import mp_tools
from astroquery.simbad import Simbad 


"""
This is the MoonPy master script! To open, you should only have to type 'import moonpy'

This package is designed to do the following:

1) download light curves from Kepler, K2 or TESS (kplr, tess)
2) generate moon models based on user specifications
3) detrend data using CoFiAM or untrendy
4) fit a model to the data using MultiNest or emcee
5) visualize the results

NOTES to Alex:
-a light curve should be an OBJECT, which you can manipulate with methods such as 'detrend', and 'fit'.
-build the skeleton first, and then start filling out the functionality.


"""


class MoonPy_LC(object):
	### this is the light curve object. You can detrend this light curve, or fit it.
	### when you initialize it, you'll either give it the times, fluxes, and errors, OR
	### you'll provide a targetID and telescope, which will allow you to download the dataset!

	def __init__(self, lc_times=None, lc_fluxes=None, lc_errors=None, targetID=None, target_type='koi', quarters='all', lc_format='pdc', coord_format='degrees', search_radius=5, sc=False, RA=None, Dec=None, telescope=None, ffi='y', lc_meta=None):
		if (lc_times != None) and (lc_fluxes != None) and (lc_errors != None):
			### implies you've supplied times, fluxes, and errorsm, so these attributes are meaningless.
			self.target = None
			self.telescope = None 
			self.meta = None


		### HANDLING FOR DOWNLOADING A LIGHT CURVE.
		elif (targetID != None) and (telescope != None):
			### implies you've selected a target you want to download.
			if telescope == 'kepler':
				### download the light curve with kplr
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, type=target_type, quarters=quarters, lc_format=lc_format, sc=sc)
			elif telescope == 'tess':
				if ffi == 'y':
					lc_times, lc_fluxes, lc_errors = eleanor_target_download(targetID)
				elif ffi == 'n':
					lc_times, lc_fluxes, lc_errors = tess_target_download(targetID)

			self.telescope = telescope
			if target_type == 'koi':
				target_name = "KOI-"+str(targetID)
			elif target_type == 'planet':
				target_name = "Kepler-"+str(targetID)
			elif target_type == 'kic':
				target_name = "KIC "+str(targetID)

			self.target = target_name 
			self.RA = Simbad.query_object(target_name)[0]['RA']
			self.Dec = Simbad.query_object(target_name)[0]['DEC']


		### HANDLING FOR A USER-SUPPLIED LIGHT CURVE.
		elif (targetID == None) and (RA != None) and (Dec != None) and (telescope != None): 
			### implies you have coordinates but not a valid target.
			if telescope == 'kepler':
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
			elif telescope == 'tess':
				lc_times, lc_fluxes, lc_errors = eleanor_coord_download(RA, Dec)



		### YOU HAVEN'T DOWNLOADED A LIGHT CURVE OR SUPPLIED ONE, SO WHAT ARE YOU DOING?
		else:
			raise Exception("You have supplied inconsistent inputs. Must be 1) lc_times, \
				lc_fluxes, and lc_errors, 2) targetID and telescope, or 3) RA, Dec and telescope.")

		self.times = lc_times 
		self.fluxes = lc_fluxes
		self.errors = lc_errors
		self.fluxes_detrend = lc_fluxes ### this will be updated when you run detrend
		self.errors_detrend = lc_errors ### this will be updated when you run detrend




	### DETRENDING!

	def detrend(self, detrend_algorithm='cofiam'):
		### optional values for method are "cofiam", "untrendy", and "george"
		if detrend_algorithm == 'cofiam':
			self.fluxes_detrend, self.errors_detrend = cofiam_detrend(self.times, self.fluxes, self.errors)
		elif detrend_algorithm == 'untrendy':
			### needs to change
			self.fluxes_detrend, self.errors_detrend = untrendy_detrend(self.times, self.fluxes, self.errors)		
		elif detrend_algorithm == 'george':
			### needs to change
			self.fluxes_detrend, self.errors_detrend = george_detrend(self.times, self.fluxes, self.errors)



	### FITTING!

	def fit(self, fitter='multinest', params=None, cost_function=None):
		### optional values for code are "multinest" and "emcee"
		if type(params) != dict:
			raise Exception("'params' must be a dictionary, with strings as the keys and priors for the values.")

		if (self.fluxes == self.fluxes_detrend) or (self.errors == self.errors_detrend):
			warnings.warn("light curve has not been detrended. It will be performed now.")
			### alternatively, just detrend!
			self.detrend()

		if fitter == 'multinest':
			mp_multinest(params, cost_function) ### outputs to a file

		elif fitter == 'emcee':
			mp_emcee(params, cost_function) ### outputs to a file





### build a function that allows you to easily make a moon model, based on user specifications.
def moon_generator(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y'):
	run_LUNA(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y')

















