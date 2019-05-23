from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
#import pymultinest
#import emcee
#import untrendy
#import sklearn
#import kplr
#import everest
#import eleanor
#import corner
#import forecaster
#import pyluna
#import tensorflow
#import cofiam
import astropy
import warnings
from astropy.io import ascii
import time
import datetime

### The packages below interface with standard packages, but within the context
### of what you want here... maybe change this usage below!
from mp_lcfind import kplr_target_download, kplr_coord_download, eleanor_target_download, eleanor_coord_download
from mp_detrend import untrendy_detrend, cofiam_detrend, george_detrend, medfilt_detrend
from mp_fit import mp_multinest, mp_emcee
from cofiam import max_order
#from pyluna import run_LUNA
import mp_tools
import traceback
from astroquery.simbad import Simbad 
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
from pyluna import run_LUNA, prepare_files


#from matplotlib import rc

#rc('font', **{'family':'serif','serif':['computer modern roman']})
#rc('text', usetex=True)

savepath = '/Users/hal9000/Documents/Software/MoonPy/saved_lcs'

"""
This is the MoonPy master script! To open, you should only have to type 'import moonpy'
This package is designed to do the following:
1) download light curves from Kepler, K2 or TESS (kplr, k2plr, tess)
2) generate moon models based on user specifications
3) detrend data using CoFiAM or untrendy
4) fit a model to the data using MultiNest or emcee
5) visualize the results
NOTES to Alex:
-a light curve should be an OBJECT, which you can manipulate with methods such as 'detrend', and 'fit'.
-build the skeleton first, and then start filling out the functionality.
"""

### STANDARD LUNA FIT PRIOR DICTIONARY. YOU MAY CHANGE INDIVIDUAL ASPECTS WHEN CALLING THE FUNCTION!
### have default param_labels, prior forms, and limits set!
### tau0 is planet specific. This will also allow you to generate the moon outside pymultinest.
param_uber_dict = {}
#param_uber_dict['tau0'] = ['uniform', (self.tau0-0.1, self.tau0+0.1)]
param_uber_dict['Rstar'] = ['loguniform', (1e6, 1e10)] ### meters 
param_uber_dict['Mstar'] = ['loguniform', (1e29, 1e33)] ### kg
param_uber_dict['q1'] = ['uniform', (0,1)]
param_uber_dict['q2'] = ['uniform', (0,1)]
param_uber_dict['Rplan'] = ['loguniform', (1e6, 1e8)]
param_uber_dict['Mplan'] = ['loguniform', (1e22, 1e30)]
param_uber_dict['bplan'] = ['uniform', (0,1)]
#param_uber_dict['Pplan'] = ['uniform', (1,3000)]
param_uber_dict['Rsat'] = ['loguniform', (1e5, 1e7)]
param_uber_dict['Msat'] = ['loguniform', (1e21, 1e27)]
param_uber_dict['sat_sma'] = ['loguniform', (1,1e3)] #### units of Rp!
param_uber_dict['sat_inc'] = ['uniform', (0,2*np.pi)]
param_uber_dict['sat_phase'] = ['uniform', (0,2*np.pi)]
param_uber_dict['sat_omega'] = ['uniform', (0,2*np.pi)]


class MoonpyLC(object):
	### this is the light curve object. You can detrend this light curve, or fit it.
	### when you initialize it, you'll either give it the times, fluxes, and errors, OR
	### you'll provide a targetID and telescope, which will allow you to download the dataset!

	def __init__(self, targetID=None, lc_times=None, lc_fluxes=None, lc_errors=None, target_type=None, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', sc=False, ffi='y', lc_meta=None, save_lc='y', tau0=None, Pplan=None):
	
		if telescope == None: # if user hasn't specified, figure it out!
			if (str(targetID).startswith("KIC")) or (str(targetID).startswith('Kepler')) or (str(targetID).startswith('kepler')) or (str(targetID).startswith('KOI')):
				telescope='kepler'
			elif str(targetID).startswith('TIC') or str(targetID).startswith('TOI'):
				telescope='tess'
			elif (str(targetID).startswith("K2")) or (str(targetID).startswith('k2')) or (str(targetID).startswith('EPIC')):
				telescope='k2'
			else:
				telescope = input('Please specify the telescope: ')

			self.telescope = telescope

		### strip off prefixes from the targetID
		if str(targetID).startswith('KIC'):
			targetID = targetID[3:]
		elif str(targetID).startswith('TIC'):
			targetID = targetID[3:]
		elif str(targetID).startswith('Kepler-'):
			targetID = targetID[7:]
		elif str(targetID).startswith('kepler-'):
			targetID = targetID[7:]
		elif str(targetID).startswith('KOI-'):
			targetID = targetID[4:]
		elif str(targetID).startswith('TOI-'):
			targetID = targetID[4:]
		elif str(targetID).startswith('K2-'):
			targetID = targetID[3:]
		elif str(targetID).startswith("EPIC"):
			targetID = targetID[4:]

		if str(targetID).startswith(' ') or str(targetID).startswith('-'):
			targetID = targetID[1:]


		### intuit whether the targetID is a 'planet' (includes a letter), a KOI (includes a decimal), or a KIC (neither).
		if target_type==None: ### not specified.
			if '.' in str(targetID) and ((telescope=='kepler') or (telescope=="Kepler")):
				target_type='koi'
			elif '.' in str(targetID) and ((telescope=='tess') or (telescope=='Tess') or (telescope=='TESS')):
				target_type='toi'
			elif (('b' in str(targetID)) or ('c' in str(targetID)) or ('d' in str(targetID)) or ('e' in str(targetID)) or ('f' in str(targetID)) or ('g' in str(targetID))) and ((telescope=='kepler') or (telescope=='Kepler')):
				target_type='planet'
			elif (telescope=='k2') or (telescope=='K2'):
				target_type='planet'
			else:
				if ((telescope == 'kepler') or (telescope=='Kepler')):
					target_type='kic'
				elif ((telescope == 'tess') or (telescope=='Tess') or (telescope=='TESS')):
					target_type = 'tic'

		print("targetID = ", targetID)
		print('target_type = ', target_type)
		print('telescope = ', telescope)

		if (type(lc_times) != type(None)) and (type(lc_fluxes) != type(None)) and (type(lc_errors) != type(None)):
			### implies you've supplied times, fluxes, and errorsm, so these attributes are meaningless.
			self.target = None
			self.telescope = None 
			self.meta = None

		### HANDLING FOR DOWNLOADING A LIGHT CURVE.
		elif (targetID != None) and (telescope != None):
			### implies you've selected a target you want to download.
			if (telescope == 'kepler') or (telescope=="Kepler") or (telescope=='k2') or (telescope=='K2'):
				### download the light curve with kplr
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, type=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
			elif (telescope == 'tess') or (telescope == "Tess") or (telescope == "TESS"):
				if ffi == 'y':
					lc_times, lc_fluxes, lc_errors = eleanor_target_download(targetID, lc_format=lc_format, sectors=quarters)
				elif ffi == 'n':
					lc_times, lc_fluxes, lc_errors = tess_target_download(targetID)

			self.telescope = telescope
			if target_type == 'koi':
				target_name = "KOI-"+str(targetID)
			elif target_type == 'toi':
				target_name = 'TOI-'+str(targetID)
			elif target_type == 'planet':
				if (telescope == 'kepler') or (telescope=="Kepler"):
					target_name = "Kepler-"+str(targetID)
				elif (telescope == 'k2') or (telescope=='K2'):
					target_name = 'K2-'+str(targetID)
			elif target_type == 'kic':
				target_name = "KIC "+str(targetID)
			elif target_type == 'tic':
				target_name = "TIC "+str(targetID)

			self.target = target_name 
			self.RA = Simbad.query_object(target_name)[0]['RA']
			self.Dec = Simbad.query_object(target_name)[0]['DEC']


		### HANDLING FOR A USER-SUPPLIED LIGHT CURVE.
		elif (targetID == None) and (RA != None) and (Dec != None) and (telescope != None): 
			### implies you have coordinates but not a valid target.
			if (telescope == 'kepler') or (telescope=='Kepler') or (telescope=='K2') or (telescope == 'k2'):
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
			elif telescope == 'tess':
				lc_times, lc_fluxes, lc_errors = eleanor_coord_download(RA, Dec)

			self.target = target_name


		### YOU HAVEN'T DOWNLOADED A LIGHT CURVE OR SUPPLIED ONE, SO WHAT ARE YOU DOING?
		else:
			raise Exception("You have supplied inconsistent inputs. Must be 1) lc_times, \
				lc_fluxes, and lc_errors, 2) targetID and telescope, or 3) RA, Dec and telescope.")


		### MAKE THEM INTO ARRAYS
		lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags)
		
		for qidx in np.arange(0,lc_times.shape[0],1):

			### remove nans
			nan_idxs = np.where(np.isfinite(lc_fluxes[qidx]) == False)[0]
			lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], nan_idxs), np.delete(lc_fluxes[qidx], nan_idxs), np.delete(lc_errors[qidx], nan_idxs), np.delete(lc_fluxes_detrend[qidx], nan_idxs), np.delete(lc_errors_detrend[qidx], nan_idxs), np.delete(lc_flags[qidx], nan_idxs)

			assert np.all(np.isfinite(lc_fluxes[qidx]))
			assert np.all(np.isfinite(lc_errors[qidx]))

			if remove_flagged == 'y':
				flag_idxs = np.where(lc_flags[qidx] != 0)[0]
				lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], flag_idxs), np.delete(lc_fluxes[qidx], flag_idxs), np.delete(lc_errors[qidx], flag_idxs), np.delete(lc_fluxes_detrend[qidx], flag_idxs), np.delete(lc_errors_detrend[qidx], flag_idxs), np.delete(lc_flags[qidx], flag_idxs)

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
		self.flags = lc_flags
		self.quarters = lc_quarters

		if save_lc == 'y':
			### write to a file!
			lcfile = open(savepath+'/'+str(target_name)+'_lightcurve.csv', mode='w')
			lcfile.write('BKJD,fluxes,errors,flags,quarter\n')
			for qidx in np.arange(0,len(self.quarters),1):
				qtq = lc_quarters[qidx]
				qtimes, qfluxes, qerrors, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_flags[qidx]
				for qt, qf, qe, qfl in zip(qtimes, qfluxes, qerrors, qflags):
					lcfile.write(str(qt)+','+str(qf)+','+str(qe)+','+str(qfl)+','+str(qtq)+'\n')
			lcfile.close()

	### DETRENDING!

	def detrend(self, dmeth='cofiam', save_lc='y', mask_transits='y', skip_ntqs='n', kernel=None, max_degree=30):
		### optional values for method are "cofiam", "untrendy", "medfilt"
		### EACH QUARTER SHOULD BE DETRENDED INDIVIDUALLY!

		try:
			self.duration_hours ### tests whether you've called get_properties yet.
		except:
			self.get_properties()

		master_detrend, master_error_detrend = [], []

		for qidx in np.arange(0,self.times.shape[0],1):
			skip_quarter = 'n'
			print('quarter = ', self.quarters[qidx])
			dtimes, dfluxes, derrors = self.times[qidx], self.fluxes[qidx], self.errors[qidx]
			print('dtimes.shape = ', dtimes.shape)

			dtimesort = np.argsort(dtimes)
			dtimes, dfluxes, derrors = dtimes[dtimesort], dfluxes[dtimesort], derrors[dtimesort]

			for ndt,dt in enumerate(dtimes):
				try:
					timediff = dtimes[ndt+1] - dtimes[ndt]
					if timediff <= 0:
						print("WARNING: timestamps are not sorted.")
				except:
					pass

			### identify in-transit idxs, and mask them!
			if mask_transits == 'y':
				### find out which transit midtimes, if any, are in this quarter
				mask_transit_idxs = []
				quarter_transit_taus = self.taus[((self.taus > np.nanmin(dtimes)) & (self.taus < np.nanmax(dtimes)))]
				for qtt in quarter_transit_taus:
					in_transit_idxs = np.where((dtimes >= qtt - self.duration_days) & (dtimes <= qtt + self.duration_days))[0]
					mask_transit_idxs.append(in_transit_idxs)
				try:
					print("transit midtimes this quarter: ", quarter_transit_taus)
					print('min, max quarter times: ', np.nanmin(dtimes), np.nanmax(dtimes))
					mask_transit_idxs = np.concatenate((mask_transit_idxs))
					print("transit in this quarter.")
				except:
					mask_transit_idxs = np.array([])
					print('no transits in this quarter.')
					if skip_ntqs == 'y':
						### skip this quarter! there are no transits present.
						fluxes_detrend, errors_detrend = dfluxes, derrors
						skip_quarter = 'y'

			elif mask_transits == 'n':
				mask_transit_idxs = None 

			if skip_quarter == 'n':
				if dmeth == 'cofiam':
					max_degree = max_order(dtimes, self.duration_days)
					print("cofiam maximum k = "+str(max_degree))
					"""
					CoFiAM is designed to preserve short-duration (i.e. moon-like) features by
					specifying a maximum order, above which the algorithm will not fit.
					This maximum degree should be calculated based on the transit duration, 
					but in practice going above k=30 is computational impractical. 
					You can specify "max_degree" below, or calculate it using the max_order function
					within cofiam.py.
					"""
					fluxes_detrend, errors_detrend = cofiam_detrend(dtimes, dfluxes, derrors, mask_idxs=mask_transit_idxs, max_degree=max_degree)

				elif dmeth == 'untrendy':
					fluxes_detrend, errors_detrend = untrendy_detrend(dtimes, dfluxes, derrors, mask_idxs=mask_transit_idxs)		

				elif dmeth == 'george':
					fluxes_detrend, errors_detrend = george_detrend(dtimes, dfluxes, derrors, mask_idxs=mask_transit_idxs)

				elif dmeth == 'medfilt':
					fluxes_detrend, errors_detrend = medfilt_detrend(dtimes, dfluxes, derrors, size=kernel, mask_idxs=mask_transit_idxs)


			### update self -- just this quarter!
			if skip_ntqs == 'n':
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
			
			lcfile = open(savepath+'/'+str(self.target)+'_lightcurve.csv', mode='w')
			lcfile.write('BKJD,fluxes,errors,fluxes_detrended,errors_detrended,flags,quarter\n')

			for qidx in np.arange(0,len(self.quarters),1):
				qtq = self.quarters[qidx]
				qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx]
				for qt, qf, qe, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags):
					lcfile.write(str(qt)+','+str(qf)+','+str(qe)+','+str(qfd)+','+str(qed)+','+str(qfl)+','+str(qtq)+'\n')
			lcfile.close()




	### FITTING!

	def fit(self, custom_param_dict=None, fitter='multinest', model='LUNA', nlive=1000):
		### optional values for code are "multinest" and "emcee"
		#if type(params) != dict:
		#	raise Exception("'params' must be a dictionary, with strings as the keys and priors for the values.")

		try:
			print(self.fluxes_detrend)
		except:
			if (self.fluxes == self.fluxes_detrend) or (self.errors == self.errors_detrend):
				warnings.warn("light curve has not been detrended. It will be performed now.")
				### alternatively, just detrend!
				self.detrend()
			else:
				self.detrend()



		### prepare seriesP.jam and plotit.f90... only needs to be done once!
		prepare_files(np.hstack(self.times))


		### the param_uber_dict is initialized at the top of this script.
		### the only standard, object-specific parameter that must be supplied is tau0.
		param_uber_dict['tau0'] = ['uniform', (self.tau0-0.1, self.tau0+0.1)]
		param_uber_dict['Pplan'] = ['uniform', (self.period-1, self.period+1)]

		if custom_param_dict != None:
			### update the parameter dictionary values!!!
			for cpdkey in custom_param_dict.keys():
				param_uber_dict[cpdkey] = custom_param_dict[cpdkey]

		global param_labels
		global param_prior_forms
		global param_limit_tuple

		param_labels = []
		param_prior_forms = []
		param_limit_tuple = []

		for pkey in param_uber_dict.keys():
			param_labels.append(pkey)
			param_prior_forms.append(param_uber_dict[pkey][0])
			param_limit_tuple.append(param_uber_dict[pkey][1])


		#param_labels = np.array(param_labels)
		#param_prior_forms = np.array(param_prior_forms)
		#param_limit_tuple = np.array(param_limit_tuple)

		print('moonpy param_labels = ', param_labels)
		print(' ')
		print('moonpy param_prior_forms = ', param_prior_forms)
		print(' ')
		print('moonpy param_limit_tuple = ', param_limit_tuple)

		#global param_labels
		#global param_prior_forms
		#global param_limit_tuple

		if fitter == 'multinest':
			mp_multinest(np.hstack(self.times), np.hstack(self.fluxes_detrend), np.hstack(self.errors_detrend), param_labels=param_labels, param_prior_forms=param_prior_forms, param_limit_tuple=param_limit_tuple, nlive=nlive, targetID=self.target) ### outputs to a file
			#mp_multinest(np.hstack(self.times), np.hstack(self.fluxes_detrend), np.hstack(self.errors_detrend), nlive=nlive, targetID=self.target) ### outputs to a file


		elif fitter == 'emcee':
			mp_emcee(params, cost_function) ### outputs to a file



	def fold(self, detrended='y'):
		### this method will phase fold your light curve. 
		fold_times = np.hstack(self.times) - np.nanmin(np.hstack(self.times)) - (0.5*self.period)
		fold_times = fold_times % self.period

		if detrended == 'y':
			plt.scatter(fold_times, np.hstack(self.fluxes_detrend), facecolors='LightCoral', edgecolors='k', s=10)
		else:
			plt.scatter(fold_times, np.hstack(self.fluxes_detrend), facecolors='LightCoral', edgecolors='k', s=10)
		plt.show()

		self.fold_times = fold_times
		self.fold_fluxes = np.hstack(self.fluxes_detrend)
		self.fold_errors = np.hstack(self.errors_detrend)




	def plot(self, facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', include_flagged='n', detrended='y', show_errors='n'):
		try:
			plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags, self.quarters

		except:
			print("WARNING: light curve has not been detrended yet!")
			plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = self.times, self.fluxes, self.errors, self.fluxes, self.errors, self.flags, self.quarters		


		### first step is to stitch the light curve together
		if type(quarters) != type('all'):
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

		plt.scatter(stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
		if show_errors == 'y':
			plt.errorbar(stitched_times, stitched_fluxes, yerr=stitched_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
		plt.xlabel('BKJD')
		plt.ylabel('Flux')
		plt.show()


	def get_properties(self):
		### first check to see if this file was already downloaded today. If not, download it!
		try:
			if (self.telescope == 'kepler') or (self.telescope == "Kepler"):
				filecreated_time = os.path.getctime('cumkois.txt')
			elif (self.telescope == 'k2') or (self.telescope == 'K2'):
				filecreated_time = os.path.getctime('cumk2ois.txt')

			current_time = time.time()
			if (current_time - filecreated_time) > 86400: ### the file is more than a day old.
				### download a new version!
				if (self.telescope == 'kepler') or (self.telescope == "Kepler"):
					os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_period,koi_period_err1,koi_period_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_prad,koi_prad_err1,koi_prad_err2,ra,dec&order=dec&format=ascii" -O "cumkois.txt"')
				elif (self.telescope == 'k2') or (self.telescope == "K2"):
					os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=k2candidates&select=epic_name,epic_candname,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,ra,dec&order=dec&format=ascii" -O "cumk2ois.txt"')
		except:
			if (self.telescope == 'kepler') or (self.telescope=='Kepler'):
				os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_period,koi_period_err1,koi_period_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_prad,koi_prad_err1,koi_prad_err2,ra,dec&order=dec&format=ascii" -O "cumkois.txt"')
			elif (self.telescope == 'k2') or (self.telescope == 'K2'):
				os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=k2candidates&select=epic_name,epic_candname,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,ra,dec&order=dec&format=ascii" -O "cumk2ois.txt"')


		### find by KICID, KOI number of planet!
		if (self.telescope == 'kepler') or (self.telescope == "Kepler"):
			cumkoi_data = ascii.read('cumkois.txt')
		elif (self.telescope == 'k2') or (self.telescope == 'K2'):
			cumkoi_data = ascii.read('cumk2ois.txt')
		cumkoi_columns = cumkoi_data.columns

		if (self.telescope == 'Kepler') or (self.telescope == 'kepler'):
			if str(self.target).startswith('Kepler-'):
				target_letter = str(self.target)[-1]
				if ' ' in self.target: ### already in the correct format, with a space between the letter.
					NEA_targetname = self.target
				else: #### there isn't a space, so add it!
					NEA_targetname = self.target[:-1]+' '+target_letter
				rowidx = np.where(cumkoi_data['kepler_name'] == NEA_targetname)[0]

			elif str(self.target).startswith('KIC'):
				NEA_targetname = int(self.target[4:])
				rowidx = np.where(cumkoi_data['kepid'] == NEA_targetname)[0]

			elif str(self.target).startswith('KOI'):
				NEA_targetname = str(self.target[4:])
				if len(NEA_targetname) == 7: ### of the form 5084.01
					NEA_targetname = 'K0'+str(NEA_targetname)
				elif len(NEA_targetname) == 6: ### of the form 163.01
					NEA_targetname = 'K00'+str(NEA_targetname)
				elif len(NEA_targetname) == 5: ### of the form 23.01
					NEA_targetname = 'K000'+str(NEA_targetname)
				elif len(NEA_targetname) == 4: ### of the form 1.01
					NEA_targetname = 'K0000'+str(NEA_targetname)
				rowidx = np.where(cumkoi_data['kepoi_name'] == NEA_targetname)[0]

		elif (self.telescope == 'k2') or (self.telescope == 'K2'):
			if str(self.target).startswith('K2-'):
				target_letter = str(self.target[-1])
				if ' ' in self.target:
					NEA_targetname = self.target
				rowidx = np.where(cumkoi_data['pl_name'] == NEA_targetname)[0]

				print("number of rows matching this description = ", len(rowidx))


		### now with the rowidx we can access the other properties we want!
		if (self.telescope == 'Kepler') or (self.telescope == 'kepler'):
			target_period, target_period_uperr, target_period_lowerr = cumkoi_data['koi_period'][rowidx], cumkoi_data['koi_period_err1'][rowidx], cumkoi_data['koi_period_err2'][rowidx]
			target_tau0, target_tau0_uperr, target_tau0_lowerr = cumkoi_data['koi_time0bk'][rowidx], cumkoi_data['koi_time0bk_err1'][rowidx], cumkoi_data['koi_time0bk_err2'][rowidx]
			target_impact, target_impact_uperr, target_impact_lowerr = cumkoi_data['koi_impact'][rowidx], cumkoi_data['koi_impact_err1'][rowidx], cumkoi_data['koi_impact_err2'][rowidx]
			target_duration, target_duration_uperr, target_duration_lowerr = cumkoi_data['koi_duration'][rowidx], cumkoi_data['koi_duration_err1'][rowidx], cumkoi_data['koi_duration_err2'][rowidx]
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = cumkoi_data['koi_ror'][rowidx], cumkoi_data['koi_ror_err1'][rowidx], cumkoi_data['koi_ror_err2'][rowidx]
			target_rp, target_rp_uperr, target_rp_lowerr = cumkoi_data['koi_prad'][rowidx], cumkoi_data['koi_prad_err1'][rowidx], cumkoi_data['koi_prad_err2'][rowidx]

		elif (self.telescope == 'k2') or (self.telescope == "K2"):
			target_period, target_period_uperr, target_period_lowerr = np.nanmedian(cumkoi_data['pl_orbper'][rowidx]), np.nanmedian(cumkoi_data['pl_orbpererr1'][rowidx]), np.nanmedian(cumkoi_data['pl_orbpererr2'][rowidx])
			target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nanmedian(cumkoi_data['pl_tranmid'][rowidx]), np.nanmedian(cumkoi_data['pl_tranmiderr1'][rowidx]), np.nanmedian(cumkoi_data['pl_tranmiderr2'][rowidx])
			target_impact, target_impact_uperr, target_impact_lowerr = np.nanmedian(cumkoi_data['pl_imppar'][rowidx]), np.nanmedian(cumkoi_data['pl_impparerr1'][rowidx]), np.nanmedian(cumkoi_data['pl_impparerr2'][rowidx])
			target_duration, target_duration_uperr, target_duration_lowerr = np.nanmedian(cumkoi_data['pl_trandur'][rowidx]), np.nanmedian(cumkoi_data['pl_trandurerr1'][rowidx]), np.nanmedian(cumkoi_data['pl_trandurerr2'][rowidx])


		### update properties!
		self.period = float(target_period)
		self.period_err = (float(target_period_lowerr), float(target_period_uperr))
		self.tau0 = float(target_tau0)
		self.tau0_err = (float(target_tau0_lowerr), float(target_tau0_uperr))
		self.impact = float(target_impact)
		self.impact_err = (float(target_impact_lowerr), float(target_impact_uperr))
		self.duration_hours = float(target_duration)
		self.duration_hours_err = (float(target_duration_lowerr), float(target_duration_uperr))
		self.duration_days = float(target_duration)/24
		self.duration_days_err = (float(target_duration_lowerr)/24, float(target_duration_uperr)/24)
		self.rprstar = float(target_rprstar)
		self.rprstar_err = (float(target_rprstar_lowerr), float(target_rprstar_uperr))
		self.rp_rearth = float(target_rp) ### earth radii
		self.rp_rjup = float(target_rp) * (R_earth.value / R_jup.value)
		self.rp_meters = self.rp_rjup * mp_tools.eq_RJup
		self.rp_rearth_err = (float(target_rp_lowerr), float(target_rp_uperr))
		self.rstar_rsol = (float(target_rp) * (1/float(target_rprstar))) * (R_earth.value / R_sun.value)
		self.rstar_meters = self.rstar_rsol * mp_tools.eq_RSun
		self.depth = self.rprstar**2

		###	identify in-transit times
		transit_midtimes = [self.tau0]
		next_transit = transit_midtimes[-1]+self.period
		while next_transit < np.nanmax(np.concatenate((self.times))):
			transit_midtimes.append(next_transit)
			next_transit = transit_midtimes[-1]+self.period
		self.taus = np.array(transit_midtimes)



### test that you can generate a moon model with pyluna from here.
### IT WORKS!
"""
ktest = MoonpyLC(targetID='Kepler-1625b', quarters=np.array([7,13,16]))
ktest.get_properties()

ktest.detrend()

prepare_files(np.hstack(ktest.times))
print('seriesP.jam and plotit.f90 have been prepared and compiled.')


ktest_modtimes, ktest_modfluxes= run_LUNA(np.hstack(ktest.times), ktest.tau0, ktest.rstar_meters, 1e31, 0.5, 0.5, ktest.rp_meters, 1e28, ktest.impact, ktest.period, mp_tools.eq_RNep, mp_tools.MNep, 25.0, 0.0, 0.0, 0.0)



#ready_to_fit = input('Ready to fit? y/n: ')
#if ready_to_fit == 'y':
#	ktest.fit()


### calculate a loglikelihood for this random test!
loglikelihood = np.nansum(-0.5 * ((ktest_modfluxes - np.hstack(ktest.fluxes_detrend)) / np.hstack(ktest.errors_detrend))**2)
print("loglikelihood = ", loglikelihood)

plt.scatter(np.hstack(ktest.times), np.hstack(ktest.fluxes_detrend), facecolors='b', alpha=0.5, s=10)
plt.scatter(ktest_modtimes, ktest_modfluxes, facecolors='LightCoral', edgecolors='k', s=10)
plt.show()
"""


#pyluna.run_LUNA(times, tau0, Rstar, Mstar, q1, q2, Rplan, Mplan, bplan, Pplan, Rsat, Msat, sat_sma, sat_inc, sat_phase, sat_omega, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='n')


"""
### build a function that allows you to easily make a moon model, based on user specifications.
def moon_generator(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y'):
	#run_LUNA(tau0, nepochs, time_from_midtime_days, cadence_minutes, noise_ppm, star_params, plan_params, sat_params, munit='kg', runit='meters', ang_unit='radians', add_noise='y', show_plots='n', print_params='n', binned_output='y')
	run_LUNA(np.hstack(self.times), tau0=self.tau0, **param_dict, add_noise='n', show_plots='y')
"""













