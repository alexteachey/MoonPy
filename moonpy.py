from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
#import pymultinest
#import emcee
#import untrendy
#import sklearn
import kplr
#import everest
#import eleanor
#import corner
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
import traceback
from astroquery.simbad import Simbad 
#from matplotlib import rc

#rc('font', **{'family':'serif','serif':['computer modern roman']})
#rc('text', usetex=True)




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

	def __init__(self, lc_times=None, lc_fluxes=None, lc_errors=None, targetID=None, target_type='koi', quarters='all', lc_format='pdc', coord_format='degrees', search_radius=5, sc=False, RA=None, Dec=None, telescope=None, ffi='y', lc_meta=None, save_lc='y', loadfile='n'):
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
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
			elif telescope == 'tess':
				lc_times, lc_fluxes, lc_errors = eleanor_coord_download(RA, Dec)

			self.target = target_name


		### YOU HAVEN'T DOWNLOADED A LIGHT CURVE OR SUPPLIED ONE, SO WHAT ARE YOU DOING?
		else:
			if loadfile == 'n':
				raise Exception("You have supplied inconsistent inputs. Must be 1) lc_times, \
					lc_fluxes, and lc_errors, 2) targetID and telescope, or 3) RA, Dec and telescope.")

			if loadfile == 'y':
				raise Exception("This functionality not yet available.")


		### MAKE THEM INTO ARRAYS
		lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags)
		
		for qidx in np.arange(0,lc_times.shape[0],1):

			### remove nans
			nan_idxs = np.where(np.isfinite(lc_fluxes[qidx]) == False)[0]
			lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], nan_idxs), np.delete(lc_fluxes[qidx], nan_idxs), np.delete(lc_errors[qidx], nan_idxs), np.delete(lc_fluxes_detrend[qidx], nan_idxs), np.delete(lc_errors_detrend[qidx], nan_idxs), np.delete(lc_flags[qidx], nan_idxs)

			assert np.all(np.isfinite(lc_fluxes[qidx]))
			assert np.all(np.isfinite(lc_errors[qidx]))

			### sort the times here!
			timesort = np.argsort(lc_times[qidx])
			lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = lc_times[qidx][timesort], lc_fluxes[qidx][timesort], lc_errors[qidx][timesort], lc_fluxes_detrend[qidx][timesort], lc_errors_detrend[qidx][timesort], lc_flags[qidx][timesort]

			for nlct,lct in enumerate(lc_times[qidx]):
				try:
					if lc_times[qidx][nlct+1] - lc_times[qidx][nlct] < 0: 
						print("times are not strictly ascending!")
				except:
					pass 

		self.times = lc_times
		self.fluxes = lc_fluxes
		self.errors = lc_errors
		#self.fluxes_detrend = lc_fluxes
		#self.errors_detrend = lc_errors
		self.flags = lc_flags
		self.quarters = lc_quarters

		if save_lc == 'y':
			### write to a file!
			lcfile = open('saved_lcs/'+str(target_name)+'_lightcurve.csv', mode='w')
			lcfile.write('BKJD,fluxes,errors,flags,quarter\n')
			for qidx in np.arange(0,len(self.quarters),1):
				qtq = lc_quarters[qidx]
				qtimes, qfluxes, qerrors, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_flags[qidx]
				for qt, qf, qe, qfl in zip(qtimes, qfluxes, qerrors, qflags):
					lcfile.write(str(qt)+','+str(qf)+','+str(qe)+','+str(qfl)+','+str(qtq)+'\n')
			lcfile.close()

	### DETRENDING!

	def detrend(self, dmeth='cofiam', save_lc='y'):
		### optional values for method are "cofiam", "untrendy", and "george"
		### EACH QUARTER SHOULD BE DETRENDED INDIVIDUALLY!

		master_detrend, master_error_detrend = [], []

		for qidx in np.arange(0,self.times.shape[0],1):
			print('quarter = ', self.quarters[qidx])
			dtimes, dfluxes, derrors = self.times[qidx], self.fluxes[qidx], self.errors[qidx]
			print('dtimes.shape = ', dtimes.shape)

			dtimesort = np.argsort(dtimes)
			dtimes, dfluxes, derrors = dtimes[dtimesort], dfluxes[dtimesort], derrors[dtimesort]

			if dmeth == 'cofiam':
				fluxes_detrend, errors_detrend = cofiam_detrend(dtimes, dfluxes, derrors)
			elif dmeth == 'untrendy':
				fluxes_detrend, errors_detrend = untrendy_detrend(dtimes, dfluxes, derrors)		
			elif dmeth == 'george':
				fluxes_detrend, errors_detrend = george_detrend(dtimes, dfluxes, derrors)

			### update self -- just this quarter!
			assert np.all(dfluxes != fluxes_detrend)
			assert np.all(derrors != errors_detrend)

			master_detrend.append(np.array(fluxes_detrend))
			master_error_detrend.append(np.array(errors_detrend))

		### this is the first initialization of the detrended fluxes.
		self.fluxes_detrend = master_detrend
		self.errors_detrend = master_error_detrend

		if save_lc == 'y':
			### overwrite the existing file!
			lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags
			
			lcfile = open('saved_lcs/'+str(self.target)+'_lightcurve.csv', mode='w')
			lcfile.write('BKJD,fluxes,errors,fluxes_detrended,errors_detrended,flags,quarter\n')

			for qidx in np.arange(0,len(self.quarters),1):
				qtq = self.quarters[qidx]
				qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx]
				for qt, qf, qe, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags):
					lcfile.write(str(qt)+','+str(qf)+','+str(qe)+','+str(qfd)+','+str(qed)+','+str(qfl)+','+str(qtq)+'\n')
			lcfile.close()




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


	def plot(self, facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', include_flagged='n', detrended='y'):
		try:
			plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags, self.quarters

		except:
			print("WARNING: light curve has not been detrended yet!")
			plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = self.times, self.fluxes, self.errors, self.fluxes, self.errors, self.flags, self.quarters		


	

		### first step is to stitch the light curve together
		if quarters != 'all':
			### means you want only selected quarters, which should be in an array!
			qstokeep_idxs = []
			for quarter in quarters:
				if quarter in self.quarters:
					qstokeep_idxs.append(np.where(quarter == self.quarters)[0][0])
			qstokeep_idxs = np.array(qstokeep_idxs)

			plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = plot_times[qstokeep_idxs], plot_fluxes[qstokeep_idxs], plot_errors[qstokeep_idxs], plot_fluxes_detrend[qstokeep_idxs], plot_errors_detrend[qstokeep_idxs], plot_flags[qstokeep_idxs], plot_quarters[qstokeep_idxs]

		if detrended == 'y':
			stitched_times, stitched_fluxes, stitched_errors, stitched_flags, stitched_quarters = np.hstack((plot_times)), np.hstack((plot_fluxes_detrend)), np.hstack((plot_errors_detrend)), np.hstack((plot_flags)), np.hstack((plot_quarters))
		else:
			stitched_times, stitched_fluxes, stitched_errors, stitched_flags, stitched_quarters = np.hstack((plot_times)), np.hstack((plot_fluxes)), np.hstack((plot_errors)), np.hstack((plot_flags)), np.hstack((plot_quarters))			

		if include_flagged=='n':
			### remove all data points with qflag != 0
			badflag_idxs = np.where(stitched_flags != 0)[0]
			stitched_times, stitched_fluxes, stitched_errors, stitched_flags = np.delete(stitched_times, badflag_idxs), np.delete(stitched_fluxes, badflag_idxs), np.delete(stitched_errors, badflag_idxs), np.delete(stitched_flags, badflag_idxs)
			assert np.all(stitched_flags == 0)

		plt.scatter(stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10)
		plt.xlabel('BKJD')
		plt.ylabel('Flux')
		plt.show()



### build a function that allows you to easily make a moon model, based on user specifications.
def moon_generator(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y'):
	run_LUNA(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y')

















