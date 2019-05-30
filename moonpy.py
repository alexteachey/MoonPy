from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import astropy
import warnings
from astropy.io import ascii
import time
import datetime
import pandas
import traceback
from astroquery.simbad import Simbad 
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
from astropy.stats import LombScargle
import os

#### BELOW ARE MOONPY PACKAGES
import mp_tools
from mp_lcfind import kplr_target_download, kplr_coord_download, eleanor_target_download, eleanor_coord_download
from mp_detrend import untrendy_detrend, cofiam_detrend, george_detrend, medfilt_detrend
from mp_batman import run_batman
from mp_fit import mp_multinest, mp_emcee
from cofiam import max_order
from pyluna import run_LUNA, prepare_files


#from matplotlib import rc

#rc('font', **{'family':'serif','serif':['computer modern roman']})
#rc('text', usetex=True)

#moonpydir = '/Users/hal9000/Documents/Software/MoonPy'
moonpydir = os.getcwd()
savepath = moonpydir+'/saved_lcs'
print("light curve savepath = ", savepath)
if os.path.exists(savepath):
	pass
else:
	print("did not exist. generating...")
	os.system('mkdir '+savepath)

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


#run_batman(all_times, RpRstar, Rstar, bplan, Pplan, tau0, q1, q2, long_peri=0, ecc=0, Mstar=None, Mplan=None, rhostar=None, rhoplan=None, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n'):
#run_LUNA(all_times, RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan, sat_sma, sat_phase, sat_inc, sat_omega, MsatMp, RsatRp, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n', **kwargs):

class MoonpyLC(object):
	### this is the light curve object. You can detrend this light curve, or fit it.
	### when you initialize it, you'll either give it the times, fluxes, and errors, OR
	### you'll provide a targetID and telescope, which will allow you to download the dataset!

	def __init__(self, targetID=None, target_type=None, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', sc=False, ffi='y', save_lc='y', load_lc='n', clobber=None):

		if (load_lc == 'n') and (clobber == None):
			### check to see if a file already exists!
			if os.path.exists(savepath+'/'+str(targetID)+'_lightcurve.csv'):
				clobber = input('light curve already exists. Clobber? y/n: ')
				if clobber == 'n':
					load_lc = 'y'
		elif (load_lc == 'n') and (clobber == 'n'):
			load_lc = 'y'
		elif (load_lc == 'n') and (clobber == 'y'):
			pass 

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

		if load_lc == 'y':
			save_lc = 'n'


		target_name = targetID 
		self.target = target_name

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
		self.telescope = telescope


		if load_lc == 'y':
			if self.target.startswith('K2') or self.target.startswith('k2'):
				self.telescope = "k2"
			elif self.target.startswith('Kepler') or self.target.startswith("kepler"):
				self.telescope = "kepler"
			else:
				telescope = input('Please specify the telescope: ')
				self.telescope = telescope 

			try:
				pandafile = pandas.read_csv(savepath+'/'+target_name+'_lightcurve.csv')
				ptimes = np.array(pandafile['BKJD'])
				pfluxes = np.array(pandafile['fluxes'])
				perrors = np.array(pandafile['errors'])
				pquarters = np.array(pandafile['quarter'])
				pflags = np.array(pandafile['flags'])
				try:
					pfluxes_detrend = np.array(pandafile['fluxes_detrended'])
					perrors_detrend = np.array(pandafile['errors_detrended'])
				except:
					pass

				unique_quarters = np.unique(pquarters)
				lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags, lc_quarters = [], [], [], [], [], [], []

				for uq in unique_quarters:
					uqidxs = np.where(pquarters == uq)[0]
					lc_times.append(ptimes[uqidxs])
					lc_fluxes.append(pfluxes[uqidxs])
					lc_errors.append(perrors[uqidxs])
					lc_flags.append(pflags[uqidxs])
					lc_quarters.append(uq)
					try:
						lc_fluxes_detrend.append(pfluxes_detrend[uqidxs])
						lc_errors_detrend.append(perrors_detrend[uqidxs])
					except:
						pass
	
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags), np.array(lc_quarters)
				self.times, self.fluxes, self.errors, self.flags, self.quarters = lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters
				try:
					lc_fluxes_detrend, lc_errors_detrend = np.array(lc_fluxes_detrend), np.array(lc_errors_detrend)
					self.fluxes_detrend, self.errors_detrend = lc_fluxes_detrend, lc_errors_detrend 
				except:
					pass

			except:
				print("could not load the light curve from file. Will download.")
				load_lc = 'n'



		### HANDLING FOR DOWNLOADING A FRESH LIGHT CURVE.
		if (load_lc=='n') and (targetID != None) and (telescope != None):
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


		### HANDLING FOR USER-SUPPLIED COORDINATES.
		elif (load_lc == 'n') and (RA != None) and (Dec != None): 
			try:	
				if (telescope == 'kepler') or (telescope=='Kepler') or (telescope=='K2') or (telescope == 'k2'):
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
				elif telescope == 'tess':
					lc_times, lc_fluxes, lc_errors = eleanor_coord_download(RA, Dec)
			except:
				telescope = input('Specify the telescope: ')
				self.telescope = telescope 
				if (telescope == 'kepler') or (telescope=='Kepler') or (telescope=='K2') or (telescope == 'k2'):
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
				elif telescope == 'tess':
					lc_times, lc_fluxes, lc_errors = eleanor_coord_download(RA, Dec)


			self.target = target_name


		### YOU HAVEN'T DOWNLOADED A LIGHT CURVE OR SUPPLIED ONE, SO WHAT ARE YOU DOING?
		elif load_lc == 'y':
			pass

		else:
			raise Exception("You have supplied inconsistent inputs. Must be 1) lc_times, \
				lc_fluxes, and lc_errors, 2) targetID and telescope, or 3) RA, Dec and telescope.")




		### BELOW THIS LINE IS HANDLING FOR EVERY DIFFERENT WAY YOU MIGHT HAVE LOADED THE LIGHT CURVE ABOVE.
		if load_lc == 'y':
			pass ### you've already turned them into arrays.
		else:
			lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags)
		
		for qidx in np.arange(0,lc_times.shape[0],1):

			### remove nans
			#if load_lc == 'n': ### if you loaded a light curve file, the NaNs will already be removed.
			nan_idxs = np.where(np.isfinite(lc_fluxes[qidx]) == False)[0]
			if len(nan_idxs) > 0:
				lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], nan_idxs), np.delete(lc_fluxes[qidx], nan_idxs), np.delete(lc_errors[qidx], nan_idxs), np.delete(lc_fluxes_detrend[qidx], nan_idxs), np.delete(lc_errors_detrend[qidx], nan_idxs), np.delete(lc_flags[qidx], nan_idxs)

			assert np.all(np.isfinite(lc_fluxes[qidx]))
			assert np.all(np.isfinite(lc_errors[qidx]))

			#if load_lc == 'n':
			if remove_flagged == 'y':
				flag_idxs = np.where(lc_flags[qidx] != 0)[0]
				if len(flag_idxs) > 0:
					lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], flag_idxs), np.delete(lc_fluxes[qidx], flag_idxs), np.delete(lc_errors[qidx], flag_idxs), np.delete(lc_fluxes_detrend[qidx], flag_idxs), np.delete(lc_errors_detrend[qidx], flag_idxs), np.delete(lc_flags[qidx], flag_idxs)

			### sort the times here!
			if load_lc == 'n': ### if you loaded a light curve, this has already been performed.
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

		### calculate which quarters have transits!
		self.find_transit_quarters()

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

	def fit(self, custom_param_dict=None, fitter='multinest', modelcode='LUNA', skip_ntqs='y', model='M', nlive=1000, nwalkers=100, nsteps=10000, resume=True, folded=False):
		### optional values for code are "multinest" and "emcee"
		#if type(params) != dict:
		#	raise Exception("'params' must be a dictionary, with strings as the keys and priors for the values.")

		### FOUR MODELS MAY BE RUN: a planet-only model (P) with no TTVs, a TTV model (T), freely fitting the transit times, 
		### a (Z) model, which gives the moon a mass but no radius, and an (M) model, which is a fully physical moon fit.

		if modelcode == "LUNA":
			folded = False ### you should not be fitting a moon model to a folded light curve!

		self.get_properties()
		self.initialize_priors(modelcode=modelcode)
		param_uber_dict = self.param_uber_dict

		if (self.fluxes == self.fluxes_detrend) or (self.errors == self.errors_detrend):
			warnings.warn("light curve has not been detrended. It will be performed now.")
			### alternatively, just detrend!
			self.detrend()


		if folded == True:
			self.fold() ### generate the folded times, in case they haven't been already.
			skip_ntqs = 'n' ### if you have a phase-folded light curve you don't need to skip non-transiting quarters!
			lc_times, lc_fluxes, lc_errors = self.fold_times, self.fold_fluxes, self.fold_errors
			### update the tau0 prior!
			self.param_uber_dict['tau0'] = ['uniform', (self.fold_tau-0.1, self.fold_tau+0.1)]

		else:
			### stacks of arrays!
			lc_times, lc_fluxes, lc_errors = self.times, self.fluxes_detrend, self.errors_detrend


		if skip_ntqs == 'y':
			### only feed in the quarters that include a transit!
			fit_times, fit_fluxes, fit_errors = [], [], []
			for qidx in np.arange(0,self.times.shape[0],1):
				qtimes, qfluxes, qerrors = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx]
				for stau in self.taus:
					if (stau >= np.nanmin(qtimes)) and (stau <+ np.nanmax(qtimes)):
						### transit in this quarter:
						fit_times.append(qtimes)
						fit_fluxes.append(qfluxes)
						fit_errors.append(qerrors)
						break
			fit_times, fit_fluxes, fit_errors = np.hstack(np.array(fit_times)), np.hstack(np.array(fit_fluxes)), np.hstack(np.array(fit_errors))

		else:
			### using all quarters!
			if folded == False:
				fit_times, fit_fluxes, fit_errors = np.hstack(lc_times), np.hstack(lc_fluxes), np.hstack(lc_errors)
			elif folded == True:
				fit_times, fit_fluxes, fit_errors = lc_times, lc_fluxes, lc_errors ### already flattened!


		### prepare seriesP.jam and plotit.f90... only needs to be done once!
		try:
			prepare_files(np.hstack(fit_times))
		except:
			prepare_files(fit_times)


		if model == 'M':
			pass
		elif (model == 'P') or (model=='T'):
			self.param_uber_dict['RsatRp'] = ['uniform', (1e-6,1e-6)]
			self.param_uber_dict['MsatMp'] = ['uniform', (1e-6,1e-6)]
			self.param_uber_dict['sat_sma'] = ['uniform', (1e-6, 1e-6)]
			self.param_uber_dict['sat_phase'] = ['uniform', (0,0)]
			self.param_uber_dict['sat_omega'] = ['uniform', (0,0)]

		if model == 'T':
			ntransits = len(self.taus)
			for i in np.arange(1,ntransits+1,1):
				taukeyname = 'tau'+str(i)
				param_uber_dict[taukeyname] = ['uniform', (self.tau0 + i*self.period -0.1, self.tau0 + i*self.period + 0.1)]

		if model == 'Z':
			param_uber_dict['RsatRp'] = ['uniform', (1e-6, 1e-6)]


		if custom_param_dict != None:
			### update the parameter dictionary values!!!
			for cpdkey in custom_param_dict.keys():
				param_uber_dict[cpdkey] = custom_param_dict[cpdkey]


		### HAS TO BE DONE HERE, SINCE YOU MAY HAVE UPDATED THE DICTIONARY!
		param_labels = []
		param_prior_forms = []
		param_limit_tuple = []

		for pkey in self.param_uber_dict.keys():
			param_labels.append(pkey)
			param_prior_forms.append(param_uber_dict[pkey][0])
			param_limit_tuple.append(param_uber_dict[pkey][1])

		self.param_labels = param_labels
		self.param_prior_forms = param_prior_forms 
		self.param_limit_tuple = param_limit_tuple


		for parlab, parprior, parlim in zip(param_labels, param_prior_forms, param_limit_tuple):
			print(parlab, parprior, parlim)


		if fitter == 'multinest':
			mp_multinest(fit_times, fit_fluxes, fit_errors, param_labels=param_labels, param_prior_forms=param_prior_forms, param_limit_tuple=param_limit_tuple, nlive=nlive, targetID=self.target, modelcode=modelcode) ### outputs to a file

		elif fitter == 'emcee':
			mp_emcee(fit_times, fit_fluxes, fit_errors, param_labels=param_labels, param_prior_forms=param_prior_forms, param_limit_tuple=param_limit_tuple, nwalkers=nwalkers, nsteps=nsteps, targetID=self.target, modelcode=modelcode, resume=resume) ### outputs to a file

		### ONE FINAL STEP -- RESTORE DEFAULT VALUES (REMOVE tau0 = folded_tau) by initializing priors again.
		self.initialize_priors(modelcode=modelcode)		



	def genLS(self, show_plot = 'y'):
		### this function generates a Lomb-Scargle Periodogram!
		LSperiods = []
		LSpowers = []
		LSfaps = []
		for qidx in np.arange(0,self.times.shape[0],1):
			qtimes, qfluxes, qerrors = self.times[qidx], self.fluxes[qidx], self.errors[qidx]
			maxperiod = 0.5 * (np.nanmax(qtimes) - np.nanmin(qtimes))
			minperiod = 0.5
			minfreq, maxfreq = 1/maxperiod, 1/minperiod
			qls = LombScargle(qtimes, qfluxes, qerrors)
			qfreq, qpower = qls.autopower(minimum_frequency=minfreq, maximum_frequency=maxfreq)
			qperiods = 1/qfreq
			qfap = qls.false_alarm_probability(qpower.max())

			if show_plot == 'y':
				plt.plot(qperiods[::-1], qpower[::-1])

			LSperiods.append(qperiods)
			LSpowers.append(qpower)
			LSfaps.append(qfap)

		if show_plot == 'y':
			plt.xscale('log')
			#plt.xlim(np.nanmin(qperiods), np.nanmax(qperiods))
			plt.xlabel('Period [days]')
			plt.ylabel('Power')
			plt.title(self.target)
			plt.show()

		LSperiods, LSpowers, LSfaps = np.array(LSperiods), np.array(LSpowers), np.array(LSfaps)

		self.LSperiods = LSperiods
		self.LSpowers = LSpowers 
		self.LSfaps = LSfaps


	def find_transit_quarters(self):
		self.get_properties()
		quarter_transit_dict = {}
		for qidx,quarter in enumerate(self.quarters):
			quarter_times = self.times[qidx]
			quarter_transit = 'n'
			for tau in self.taus:
				if (tau >= np.nanmin(quarter_times)) and (tau <= np.nanmax(quarter_times)):
					### there's a transit in this quarter!
					quarter_transit = 'y'
					break
			if quarter_transit == 'n':
				quarter_transit_dict[quarter] = 'n'
			else:
				quarter_transit_dict[quarter] = 'y'
				
		self.quarter_transit_dict = quarter_transit_dict 


	def fold(self, detrended='y'):
		### this method will phase fold your light curve. 
		### first tau in the time series:
		self.get_properties()
		first_tau = self.tau0
		ftidx = 0
		while first_tau < np.nanmin(np.hstack(self.times)):
			ftidx += 1
			first_tau = self.taus[ftidx]


		fold_times = np.hstack(self.times) - np.nanmin(np.hstack(self.times)) - (0.5*self.period)
		#fold_times = np.hstack(self.times) - first_tau - (0.5*self.period)
		fold_times = fold_times % self.period
		fold_first_tau = (first_tau - np.nanmin(np.hstack(self.times)) - (0.5*self.period)) % self.period
		print ('fold_first_tau = ', fold_first_tau)
		

		if detrended == 'y':
			fold_fluxes = np.hstack(self.fluxes_detrend)
			fold_errors = np.hstack(self.errors_detrend)
		else:
			fold_fluxes = np.hstack(self.fluxes)
			fold_errors = np.hstack(self.errors)

		#fold_sort = np.argsort(fold_times)
		#fold_fluxes = np.hstack(self.fluxes_detrend)[fold_sort]
		#fold_errors = np.hstack(self.errors_detrend)[fold_sort]
		#fold_times = fold_times[fold_sort]

		self.fold_times = fold_times
		self.fold_fluxes = fold_fluxes
		self.fold_errors = fold_errors
		self.fold_tau = fold_first_tau





	def initialize_priors(self, modelcode):
		param_uber_dict = {}

		param_uber_dict['RpRstar'] = ['uniform', (0, 1)]
		param_uber_dict['rhostar'] = ['uniform', (1, 1e6)] ## roughly the density of betelgeuse to the density of Mercury.
		param_uber_dict['bplan'] = ['uniform', (0,1)]
		param_uber_dict['q1'] = ['uniform', (0,1)]
		param_uber_dict['q2'] = ['uniform', (0,1)]
		param_uber_dict['rhoplan'] = ['uniform', (1, 1e6)]
		if modelcode == "LUNA":
			### these parameters are only relevant for LUNA!
			param_uber_dict['sat_phase'] = ['uniform', (0,2*np.pi)]
			param_uber_dict['sat_inc'] = ['uniform', (-1,3)] ### actually cos(i_s), natively.
			param_uber_dict['sat_omega'] = ['uniform', (-np.pi,np.pi)]
			param_uber_dict['MsatMp'] = ['uniform', (0, 1)]
			param_uber_dict['RsatRp'] = ['uniform', (-1, 1)]
		if modelcode == 'batman':
			### these parameters are only relevant for batman!
			param_uber_dict['Rstar'] = ['loguniform', (1e6,1e11)]
			param_uber_dict['long_peri'] = ['uniform', (0,4*np.pi)] ### longitude of periastron is the sum of the ascending node  (0-2pi) and the argument of periaspe (also 0-2pi).
			param_uber_dict['ecc'] = ['uniform', (0,1)]

		### THE FOLLOWING PARAMETERS ARE PLANET-SPECIFIC.
		self.get_properties()
		param_uber_dict['tau0'] = ['uniform', (self.tau0-0.1, self.tau0+0.1)]
		param_uber_dict['Pplan'] = ['uniform', (self.period-1, self.period+1)]
		if modelcode == "LUNA":
			param_uber_dict['sat_sma'] = ['uniform', (2,7.897*(self.period**(2/3)))]

		self.param_uber_dict = param_uber_dict

		### OK TO LEAVE HERE, SO LONG AS IT ALSO STAYS WITHIN THE FIT() FUNCTION.
		param_labels = []
		param_prior_forms = []
		param_limit_tuple = []

		for pkey in self.param_uber_dict.keys():
			param_labels.append(pkey)
			param_prior_forms.append(param_uber_dict[pkey][0])
			param_limit_tuple.append(param_uber_dict[pkey][1])

		self.param_labels = param_labels
		self.param_prior_forms = param_prior_forms 
		self.param_limit_tuple = param_limit_tuple






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





	###############################
	### VISUALIZATION FUNCTIONS ###
	###############################


	def plot_lc(self, facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', folded='n', include_flagged='n', detrended='y', show_errors='n'):
		### THIS FUNCTION PLOTS THE LIGHT CURVE OBJECT.
		try:
			plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags, self.quarters

		except:
			print("WARNING: light curve has not been detrended yet!")
			detrended = 'n'
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


		if folded == 'n':
			plt.scatter(stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			if show_errors == 'y':
				plt.errorbar(stitched_times, stitched_fluxes, yerr=stitched_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
		elif folded == 'y':
			try:
				self.fold(detrended=detrended)
			except:
				self.get_properties()
				self.fold()
			plt.scatter(self.fold_times, self.fold_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			if show_errors == 'y':
				plt.errorbar(self.fold_times, self.fold_fluxes, yerr=self.fold_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
		plt.xlabel('BKJD')
		plt.ylabel('Flux')
		try:
			plt.title(str(self.target))
		except:
			pass
		plt.show()



	def plot_corner(self, fitter='emcee', modelcode='batman', burnin_pct=0.1):
		import corner

		### THIS FUNCTION GENERATES A CORNER PLOT BASED ON YOUR MODEL FITS.
		if fitter == 'multinest':
			### use this to generate a corner plot from the fit results.
			fit_resultsdir = moonpydir+'/MultiNest_fits/'+str(self.target)+'/chains'
			PEWfile = np.genfromtxt(fit_resultsdir+'/'+str(self.target)+'post_equal_weights.dat')

			json_file = open(fit_resultsdir+'/'+str(self.target)+'_params.json', mode='r')
			json_params = json_file.readline()
			json_params = json_params.split(',')

			PEWdict = {}
			for njpar, jpar in enumerate(json_params):
				while jpar.startswith(' '):
					jpar = jpar[1:]
				while jpar[-1] == ' ':
					jpar = jpar[:-1]
				while jpar.startswith('"'):
					jpar = jpar[1:]
				while jpar[-1] == '"':
					jpar = jpar[:-1]
				PEWdict[jpar] = PEWfile.T[njpar]

			### as a test, just generate a simple histogram
			for param in PEWdict.keys():
				n, bins, edges = plt.hist(PEWdict[param], bins=50, facecolor='green', edgecolor='k', alpha=0.7)
				plt.title(param)
				plt.show()

		elif fitter == 'emcee':
			if modelcode=='batman':
				chainsdir = moonpydir+'/emcee_fits/batman/'+str(self.target)+'/chains'
			elif modelcode=='LUNA':
				chainsdir=moonpydir+'/emcee_fits/LUNA/'+str(self.target)+'/chains'
			samples = np.genfromtxt(chainsdir+'/'+str(self.target)+'_mcmc_chain.txt')
			sample_shape = samples.shape
			samples = samples[int(burnin_pct*sample_shape[0]):,1:]


		self.initialize_priors(modelcode=modelcode)
		fig = corner.corner(samples, labels=self.param_labels)
		plt.savefig(chainsdir+'/'+str(self.target)+"_corner.png")
		plt.close()



	def plot_bestmodel(self, fitter, modelcode, folded=False, burnin_pct=0.1):
		### THIS FUNCTION PLOTS YOUR BEST FIT LIGHT CURVE MODEL OVER THE DATA.
		if folded == True:
			self.fold()

		if modelcode == "LUNA":
			folded = False ### should not be generating a folded light curve for a moon fit.

		self.initialize_priors(modelcode=modelcode)

		if fitter == 'emcee':
			if modelcode=='batman':
				chainsdir = moonpydir+'/emcee_fits/batman/'+str(self.target)+'/chains'
			elif modelcode=='LUNA':
				chainsdir=moonpydir+'/emcee_fits/LUNA/'+str(self.target)+'/chains'
			samples = np.genfromtxt(chainsdir+'/'+str(self.target)+'_mcmc_chain.txt')
			sample_shape = samples.shape
			samples = samples[int(burnin_pct*sample_shape[0]):,1:]

			### 
			best_fit_dict = {}
			for npar, parlab in enumerate(self.param_labels):
				best_fit_dict[parlab] = np.nanmedian(samples.T[npar])
			print("best fit values: ")
			for parkey in best_fit_dict.keys():
				print(parkey, ' = ', best_fit_dict[parkey])

			if modelcode == 'batman':
				### use batman to generate a model!!!
				if folded == True:
					batman_times, batman_fluxes = run_batman(self.fold_times, **best_fit_dict, add_noise='n', show_plots='n')
					plt.scatter(np.hstack(self.fold_times), np.hstack(self.fluxes_detrend), facecolor='LightCoral', edgecolor='k')
				else:
					batman_times, batman_fluxes = run_batman(self.times, **best_fit_dict, add_noise='n', show_plots='n')					
					plt.scatter(np.hstack(self.times), np.hstack(self.fluxes_detrend), facecolor='LightCoral', edgecolor='k')
				batman_sort = np.argsort(batman_times)
				batman_times, batman_fluxes = batman_times[batman_sort], batman_fluxes[batman_sort]
				plt.plot(batman_times, batman_fluxes, c='g', linewidth=2)
				plt.show()


