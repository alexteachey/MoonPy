from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
#import astropy
import warnings
from astropy.io import ascii
from astropy.time import Time
import time
#import datetime
import pandas
import traceback
from astroquery.simbad import Simbad 
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
from astropy.timeseries import LombScargle
import socket 

#### BELOW ARE MOONPY PACKAGES
from mp_tools import *
from mp_lcfind import *
from mp_detrend import untrendy_detrend, cofiam_detrend, george_detrend, medfilt_detrend, polyAM_detrend
from mp_batman import run_batman
from mp_fit import mp_multinest, mp_emcee
from mp_plotter import *
from cofiam import max_order
#from pyluna import run_LUNA, prepare_files
from pyluna import prepare_files
from mp_tpf_examiner import *
from scipy.interpolate import interp1d 



plt.rcParams["font.family"] = 'serif'

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/_mp_manipulation.py')]



def fit(self, custom_param_dict=None, fitter='multinest', modelcode='LUNA', segment='y', segment_length=500, skip_ntqs='y', model='M', nlive=1000, nwalkers=100, nsteps=10000, resume=True, folded=False):
	### optional values for code are "multinest" and "emcee"

	### FOUR MODELS MAY BE RUN: a planet-only model (P) with no TTVs, a TTV model (T), freely fitting the transit times, 
	### a (Z) model, which gives the moon a mass but no radius, and an (M) model, which is a fully physical moon fit.

	if modelcode == "LUNA":
		folded = False ### you should not be fitting a moon model to a folded light curve!

	##### LIGHT CURVES NEED TO HAVE THE DATA POINTS WELL OUTSIDE THE TRANSITS REMOVED,
	###### OTHERWISE IT TAKES WAAAAAAY TOO LONG TO CHEW ON THIS. (11/25/2020)

	if self.telescope.lower() != 'user':
		self.get_properties(locate_neighbor='n')
	self.initialize_priors(modelcode=modelcode)
	param_uber_dict = self.param_uber_dict

	if self.telescope.lower() == 'user':
		detrend_option = input('WARNING: lc has not been detrended. Perform now? y/n: ')
		if detrend_option == 'y':
			self.detrend()
		else:
			self.fluxes_detrend = self.fluxes
			self.errors_detrend = self.errors 
	else:
		if (self.fluxes == self.fluxes_detrend) or (self.errors == self.errors_detrend):
			warnings.warn("light curve has not been detrended. It will be performed now.")
			### alternatively, just detrend!
			self.detrend()


	if folded == True:
		self.fold() ### generate the folded times, in case they haven't been already.
		skip_ntqs = 'n' ### if you have a phase-folded light curve you don't need to skip non-transiting quarters!
		lc_times, lc_fluxes, lc_errors = self.fold_times, self.fold_fluxes, self.fold_errors
		### update the tau0 prior!
		self.param_uber_dict['tau0'] = ['uniform', (self.fold_tau-0.1, self.fold_tau+0.1)] ### establishing the prior.

	else:
		### stacks of arrays!
		lc_times, lc_fluxes, lc_errors = self.times, self.fluxes_detrend, self.errors_detrend


	if skip_ntqs == 'y':
		### only feed in the quarters that include a transit!
		if len(self.quarters) == 1:
			fit_times, fit_fluxes, fit_errors = self.times, self.fluxes_detrend, self.errors_detrend
		else:
			fit_times, fit_fluxes, fit_errors = [], [], []
			for qidx in np.arange(0,self.times.shape[0],1):
				qtimes, qfluxes, qerrors = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx]
				for stau in self.taus:
					if (stau >= np.nanmin(qtimes)) and (stau <+ np.nanmax(qtimes)):
						### transit in this quarter:

						if segment == 'y':
							#### we're going to segment the light curve, by including only n=segment_length data points on either side of tau
							##### THIS WILL GREATLY REDUCE THE COMPUTATION TIME
							stau_idx = np.nanargmin(np.abs(stau - qtimes)) #### closest index within qtimes for the transit time.
							segment_idxs = np.arange(int(stau_idx - np.ceil(segment_length/2)), int(stau_idx + np.ceil(segment_length/2)), 1)
							qtimes, qfluxes, qerrors = qtimes[segment_idxs], qfluxes[segment_idxs], qerrors[segment_idxs]

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







	if model == 'M':
		pass

	elif (model == 'P') or (model=='T'):
		### THESE ARE INPUTS TO THE MODEL BUT SHOULD NOT BE FREE PARAMETERS!
		self.param_uber_dict['RsatRp'] = ['fixed', 1e-8]
		self.param_uber_dict['MsatMp'] = ['fixed', 1e-8]
		self.param_uber_dict['sat_sma'] = ['fixed', 1000]
		self.param_uber_dict['sat_inc'] = ['fixed', 0]
		self.param_uber_dict['sat_phase'] = ['fixed', 0]
		self.param_uber_dict['sat_omega'] = ['fixed', 0]

	if model == 'T':
		transitnum = 0
		for i in np.arange(0, len(self.taus),1):
			#### verify that the tau is within your quarters!
			for qidx,q in enumerate(self.quarters):
				if self.taus[i] >= np.nanmin(self.times[qidx]) and self.taus[i] <= np.nanmax(self.times[qidx]):
					### means the transit is in your data and can be fit.
					taukeyname = 'tau'+str(transitnum)
					transitnum += 1
					param_uber_dict[taukeyname] = ['uniform', (self.taus[i] - 0.1, self.taus[i] + 0.1)]

	if model == 'Z':
		param_uber_dict['RsatRp'] = ['fixed', 1e-8]


	### prepare seriesP.jam and plotit.f90... only needs to be done once!
	if model == "M" or model == "P" or model == "Z":
		ntaus = 1
	elif model == 'T':
		ntaus = 0
		for paramkey in param_uber_dict.keys():
			if paramkey.startswith('tau'):
				ntaus += 1

	if model == "M":
		nparamorig = 14
		nparam = 14
		nvars = 14 ### fitting all the parameters!
	elif model == "P":
		nparamorig = 14 ### all these inputs must still be present, even if some of them are fixed at zero!
		nparam = 14
		nvars = 8  ### RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan
	elif model == 'T':
		nparamorig = 14
		nparam = nparamorig + (ntaus-1) ### tau0 is a STANDARD nparamorig input... every additional tau needs to be counted.
		nvars = 8 + (ntaus-1) ### standard P model variables plus all additional taus.
	elif model == 'Z':
		nparam = 14
		nparamorig = 14
		nvars = 13 ### not fitting Rsat/Rp 

	try:
		prepare_files(np.hstack(fit_times), ntaus, nparam, nparamorig)
	except:
		prepare_files(fit_times, ntaus, nparam, nparamorig)


	if custom_param_dict != None:
		### update the parameter dictionary values!!!
		for cpdkey in custom_param_dict.keys():
			param_uber_dict[cpdkey] = custom_param_dict[cpdkey]


	### HAS TO BE DONE HERE, SINCE YOU MAY HAVE UPDATED THE DICTIONARY!
	param_labels = []
	param_prior_forms = []
	param_limit_tuple = []

	for pkey in self.param_uber_dict.keys():
		param_labels.append(pkey) ### name of the prior 
		param_prior_forms.append(param_uber_dict[pkey][0]) ### the shape of the prior
		param_limit_tuple.append(param_uber_dict[pkey][1]) ### limits of the prior

	self.param_labels = param_labels
	self.param_prior_forms = param_prior_forms 
	self.param_limit_tuple = param_limit_tuple


	for parlab, parprior, parlim in zip(param_labels, param_prior_forms, param_limit_tuple):
		print(parlab, parprior, parlim)


	if fitter == 'multinest':
		mp_multinest(fit_times, fit_fluxes, fit_errors, param_dict=self.param_uber_dict, nlive=nlive, targetID=self.target, modelcode=modelcode, model=model, nparams=nvars)

	elif fitter == 'emcee':
		mp_emcee(fit_times, fit_fluxes, fit_errors, param_dict=self.param_uber_dict, nwalkers=nwalkers, nsteps=nsteps, targetID=self.target, modelcode=modelcode, model=model, resume=resume, nparams=nvars) ### outputs to a file

	### ONE FINAL STEP -- RESTORE DEFAULT VALUES (REMOVE tau0 = folded_tau) by initializing priors again.
	self.initialize_priors(modelcode=modelcode)		






def prep_for_CNN(self, save_lc='y', window=6, cnn_len=493, exclude_neighbors='y', flag_neighbors='y', show_plot='n', extra_path_info=None, cnnlc_path=moonpydir+'/cnn_lcs'):
	### this function will produce an arrayy that's ready to be fed to a CNN for moon finding. 
	if os.path.exists(cnnlc_path) == False:
		os.system('mkdir '+cnnlc_path)

	localpath_list = []	

	for taunum,tau in enumerate(self.taus):
		cnn_min, cnn_max = tau-window, tau+window
		print('cnn_min, cnn_max = ', cnn_min, cnn_max)
		cnn_idxs = np.where((np.hstack(self.times) >= cnn_min) & (np.hstack(self.times) <= cnn_max))[0]
		### grab the times, fluxes, and errors
		cnn_times, cnn_fluxes, cnn_errors, cnn_fluxes_detrend, cnn_errors_detrend = np.hstack(self.times)[cnn_idxs], np.hstack(self.fluxes)[cnn_idxs], np.hstack(self.errors)[cnn_idxs], np.hstack(self.fluxes_detrend)[cnn_idxs], np.hstack(self.errors_detrend)[cnn_idxs]
		if len(cnn_times) < cnn_len: ### the size of the array you need to feed into the CNN.
			### not enough times 
			continue
		while len(cnn_times) > cnn_len:
			### shave off from the front end.
			cnn_times, cnn_fluxes, cnn_errors, cnn_fluxes_detrend, cnn_errors_detrend = cnn_times[1:], cnn_fluxes[1:], cnn_errors[1:], cnn_fluxes_detrend[1:], cnn_errors_detrend[1:]
			if len(cnn_times) > cnn_len:
				### shave off at the back end.
				cnn_times, cnn_fluxes, cnn_errors, cnn_fluxes_detrend, cnn_errors_detrend = cnn_times[:-1], cnn_fluxes[:-1], cnn_errors[:-1], cnn_fluxes_detrend[:-1], cnn_errors_detrend[:-1]
		assert len(cnn_times) == cnn_len

		if exclude_neighbors == 'y':
			neighbor_contam = 'n'

			for neighbor in self.neighbor_dict.keys():
				neighbor_taus = self.neighbor_dict[neighbor].taus
				if np.any((neighbor_taus >= np.nanmin(cnn_times)) & (neighbor_taus <= np.nanmax(cnn_times))):
					### there's another transiting planet in your window!
					neighbor_contam = 'y'
			if neighbor_contam == 'y':
				continue

		if flag_neighbors=='y':
			flag_idxs = []
			for neighbor in self.neighbor_dict.keys():
				neighbor_taus = self.neighbor_dict[neighbor].taus
				for nt in neighbor_taus:
					if (nt >= np.nanmin(cnn_times)) and (nt <= np.nanmax(cnn_times)):
						ntidxs = np.where((cnn_times >= nt-(2.5*self.neighbor_dict[neighbor].duration_days)) & (cnn_times <= nt+(2.5*self.neighbor_dict[neighbor].duration_days)))[0]
						flag_idxs.append(ntidxs)
			#try:
			if len(flag_idxs) > 0:
				flag_idxs = np.unique(np.hstack(np.array(flag_idxs)))
			else:
				flag_idxs = np.array([])

			flag_array = np.zeros(shape=len(cnn_times))
			try:
				flag_array[flag_idxs] = 1
			except:
				pass
		### at this point you've excluded neighbors (if selected; by default) and made the light curve the right size. Save it!
		if flag_neighbors == 'y':
			cnn_stack = np.vstack((cnn_times, cnn_fluxes, cnn_errors, cnn_fluxes_detrend, cnn_errors_detrend, flag_array))
		else:
			cnn_stack = np.vstack((cnn_times, cnn_fluxes, cnn_errors, cnn_fluxes_detrend, cnn_errors_detrend))

		if type(extra_path_info) == type(None):
			localpath = cnnlc_path+'/'+self.target+"_transit"+str(taunum)+'_cnnstack.npy'
		elif type(extra_path_info) != type(None):
			localpath = cnnlc_path+'/'+self.target+'_transit'+str(taunum)+'_'+extra_path_info+'_cnnstack.npy'
		localpath_list.append(localpath)

		np.save(localpath, cnn_stack)

		if show_plot == 'y':
			plt.scatter(cnn_stack[0], cnn_stack[3], facecolor='LightCoral', edgecolor='k', s=10)
			if flag_neighbors == 'y':
				neighbor_transit_idxs = np.where(flag_array == 1)[0]
				plt.scatter(cnn_stack[0][neighbor_transit_idxs], cnn_stack[3][neighbor_transit_idxs], c='g', marker='x', s=15)
			if (self.telescope.lower() == 'kepler') or (self.telescope.lower() == 'k2'):
				plt.xlabel("BKJD")
			elif self.telescope.lower() == 'tess':
				plt.xlabel("BTJD")
			plt.ylabel('Flux')
			plt.show()

	return localpath_list 



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
	self.get_properties(locate_neighbor='n')
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



### DETRENDING!

def detrend(self, dmeth='cofiam', save_lc='y', mask_transits='y', mask_neighbors='y', skip_ntqs='n', medfilt_kernel_transit_multiple=5, GP_kernel='ExpSquaredKernel', GP_metric=1.0, max_degree=30, use_mazeh='y'):
	exceptions_raised = 'n'

	self.dmeth=dmeth
	### make sure you get_neighbors() first!
	self.get_neighbors()

	if self.telescope.lower() != 'kepler':
		use_mazeh == 'n' ### mazeh is only for Kepler targets!

	if mask_neighbors == 'y':
		mask_transits = 'y' 
	### optional values for detrending method (dmeth) are "cofiam", "untrendy", "medfilt", "george", and "polyAM"
	### EACH QUARTER SHOULD BE DETRENDED INDIVIDUALLY!

	master_detrend_model, master_detrend, master_error_detrend, master_flags_detrend = [], [], [], []

	nquarters = len(self.quarters)

	all_quarter_mask_transit_idxs = []			
	for qidx in np.arange(0,nquarters,1):
		skip_quarter = 'n'
		print('quarter = ', self.quarters[qidx])
		if nquarters != 1:
			dtimes, dfluxes, derrors, dflags = self.times[qidx], self.fluxes[qidx], self.errors[qidx], self.flags[qidx]
		elif nquarters == 1:
			dtimes, dfluxes, derrors, dflags = self.times, self.fluxes, self.errors, self.flags
		
		if dtimes.shape[0] == 0:
			exceptions_raised = 'y'
			fluxes_detrend, errors_detrend = dfluxes, derrors
			flags_detrend = np.linspace(2097152,2097152,len(fluxes_detrend))
			master_detrend.append(np.array(fluxes_detrend))
			master_error_detrend.append(np.array(errors_detrend))
			master_flags_detrend.append(np.array(flags_detrend))
			continue

		dtimesort = np.argsort(dtimes)
		dtimes, dfluxes, derrors, dflags = dtimes[dtimesort], dfluxes[dtimesort], derrors[dtimesort], dflags[dtimesort]

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

			if use_mazeh == 'y':

				print("Using TTV Catalog to identify transit times for detrending...")
				### taus should be calculated based on the Mazeh table.
				mazeh = pandas.read_csv(moonpydir+'/Table3_O-C.csv')
				maz_koi = np.array(mazeh['KOI']).astype(str)
				maz_epoch = np.array(mazeh['n'])
				maz_reftime = np.array(mazeh['ref_time'])
				maz_reftime = maz_reftime + 2454900 - 2454833 #### adjusts into BKJD -- a difference of 67 days from the Holczer catalog.	
				maz_OCmin = np.array(mazeh['O-C_min'])

				### find the indices of the target!
				target_idxs = np.where(self.target[4:] == maz_koi)[0]
				if len(target_idxs) == 0:
					print('no TTV record under the name '+str(self.target)+'. Checking aliases...')
					### check for aliases -- your self.target may not be a KOI natively!
					try:
						for planet_alias in self.aliases:
							if planet_alias.lower().startswith('koi'):
								target_idxs = np.where(planet_alias[4:] == maz_koi)[0]
								print('found a TTV record for alias: ', planet_alias)
								break
					except:
						print('Could not find any TTV record. Computing transit times from Tau0 and orbital period.')
						pass

				if len(target_idxs) == 0:
					self.mask_taus = self.taus ### it's not in the catalog. Keep moving!

				else:
					target_epochs = [] ### a list of all the epochs in the Mazeh catalog.
					target_reftimes = [] ### ditto for reftimes
					target_OCmins = [] #### ditto for O - C (minutes)
					target_OCdays = [] #### converting above to days
					for tidx in target_idxs: ### all the indices for this planet -- have to loop through to grab the relevant values.
						try:
							tepoch = int(maz_epoch[tidx])
							treftime = float(maz_reftime[tidx])
							tOCmin = float(maz_OCmin[tidx])
							tOCday = tOCmin / (60 * 24)
						except:
							traceback.print_exc()

						target_epochs.append(tepoch)
						target_reftimes.append(treftime)
						target_OCmins.append(tOCmin)
						target_OCdays.append(tOCday)

					target_epochs, target_reftimes, target_OCmins, target_OCdays = np.array(target_epochs), np.array(target_reftimes), np.array(target_OCmins), np.array(target_OCdays)

					#print('target_reftimes = ', target_reftimes)

					### interpolate OCmins! to grab missing transit times -- have to do this in case Mazeh has not cataloged them all.
					all_target_epochs = np.arange(np.nanmin(target_epochs), np.nanmax(target_epochs)+1, 1)
					all_target_reftimes = [] ### a list of all the transit times (basically linear ephemeris, from which the offset is calculated)
					all_target_OCmins = []
					all_target_OCdays = []

					OCmin_interp = interp1d(target_epochs, target_OCmins, kind='quadratic', bounds_error=False, fill_value='extrapolate')
					reftime_interp = interp1d(target_epochs, target_reftimes, kind='linear', bounds_error=False, fill_value='extrapolate')

					### the number of epochs (integer transit numbers) in "all_target_epochs" will *always* be greater or equal to number in target_epochs
					#### target epochs is the # of epochs recorded in Mazeh, all_target_epochs is every possible one from first transit to last.
					##### it appears that Mazeh DOES skip numbers when there is a missing transit, so this *ought* to work.

					for ate in all_target_epochs:
						if ate in target_epochs:
							#### this value is already recorded by mazeh
							epoch_idx = np.where(target_epochs == ate)[0]
							epoch_reftime = target_reftimes[epoch_idx]
							epoch_OCmins = target_OCmins[epoch_idx]
							epoch_OCdays = target_OCdays[epoch_idx]

						elif ate not in target_epochs:
							#### you still want it! so you're going to interpolate it.
							epoch_reftime = reftime_interp(ate)
							epoch_OCmins = OCmin_interp(ate)
							epoch_OCdays = epoch_OCmins / (60 * 24)

						all_target_reftimes.append(epoch_reftime)
						all_target_OCmins.append(epoch_OCmins)
						all_target_OCdays.append(epoch_OCdays)

					all_target_reftimes, all_target_OCmins, all_target_OCdays = np.array(all_target_reftimes), np.array(all_target_OCmins), np.array(all_target_OCdays)
					assert len(all_target_epochs) == len(all_target_reftimes)	
					assert len(all_target_epochs) >= len(target_epochs)

					### FINALLY, calculate the taus!
					mazeh_taus = all_target_reftimes + all_target_OCdays 
					self.mask_taus = mazeh_taus
					if len(self.mask_taus.shape) > 1:
						self.mask_taus = self.mask_taus.reshape(self.mask_taus.shape[0])

			elif use_mazeh == 'n':
				self.mask_taus = self.taus 

			quarter_transit_taus = self.mask_taus[((self.mask_taus > np.nanmin(dtimes)) & (self.mask_taus < np.nanmax(dtimes)))]

			for qtt in quarter_transit_taus:
				#### FOR EACH TRANSIT IN THIS QUARTER.... 

				in_transit_idxs = np.where((dtimes >= qtt - 2.5*self.duration_days) & (dtimes <= qtt + 2.5*self.duration_days))[0]
				mask_transit_idxs.append(in_transit_idxs)
			
			try:
				mask_transit_idxs = np.concatenate(mask_transit_idxs)
			except:
				print('mask_transit_idxs could not be concatenated (probably not needed).')


			#print("BEFORE: ")
			#print('mask_transit_idxs = ', mask_transit_idxs)
			#print("mask transit times = ", dtimes[mask_transit_idxs])
			mask_transit_idxs = mask_transit_idxs

			### add neighbor transit times to mask_transit_idxs.
			try:
				sntt = self.neighbor_transit_times # don't need to print them, just see if they're there!
			except:
				try:
					self.mystery_solver(self.tau0, self.period, self.duration_hours)
				except:
					pass

			if (mask_neighbors == 'y') and (len(self.neighbors) > 0):

				neighbor_transit_idxs = np.where( (self.neighbor_transit_times > np.nanmin(dtimes)) & (self.neighbor_transit_times < np.nanmax(dtimes)) )[0]	

				if len(neighbor_transit_idxs) > 0:
					print("NEIGHBORS IN THIS SEGMENT!")	
					#### THERE WON'T BE THAT MANY, JUST USE A FOR LOOP TO AVOID WEIRD FORMATTING ISSUES
					try:
						mask_transit_idxs = mask_transit_idxs.tolist()
					except:
						pass
					n_neighbor_idxs = 0
					for nti in neighbor_transit_idxs:
						mask_transit_idxs.append(nti)	
						n_neighbor_idxs += 1
					print('appended '+str(n_neighbor_idxs)+' neighbor transit data points for masking.')	

			#try:
			print('min, max quarter times: ', np.nanmin(dtimes), np.nanmax(dtimes))
			
			if len(quarter_transit_taus) > 0:
				print(str(len(quarter_transit_taus))+" transit(s) in this quarter.")

			if len(quarter_transit_taus) == 1:
				try:
					mask_transit_idxs = np.concatenate(mask_transit_idxs)
				except:
					print('Concatenate not necessary')
				mask_transit_idxs = np.array(mask_transit_idxs)				

			if len(quarter_transit_taus) > 1:
					try:
						mask_transit_idxs = np.concatenate(mask_transit_idxs)
					except:
						print('Concatenate not necessary')
					mask_transit_idxs = np.array(mask_transit_idxs)

			elif len(quarter_transit_taus) < 1:

				mask_transit_idxs = np.array([])
				print('no transits in this quarter.')
				if skip_ntqs == 'y':
					### skip this quarter! there are no transits present.
					fluxes_detrend, errors_detrend, flags_detrend = dfluxes, derrors, dflags
					skip_quarter = 'y'

			try:
				np.concatenate(mask_transit_idxs)
			except:
				pass
			mask_transit_idxs = np.unique(mask_transit_idxs)


			#### remove out of bounds transit idxs:
			OoB_mask_transit_idxs = np.where(mask_transit_idxs >= len(dtimes))[0]
			mask_transit_idxs = np.delete(mask_transit_idxs, OoB_mask_transit_idxs)
			

		elif mask_transits == 'n':
			mask_transit_idxs = np.array([])

		print(' ')
		#print("AFTER: ")
		#print('mask_transit_idxs = ', mask_transit_idxs)
		if len(mask_transit_idxs) > 0:
			print('len(mask_transit_idxs) = '+str(len(mask_transit_idxs)))
			#print("mask transit times = ", dtimes[mask_transit_idxs])			


		all_quarter_mask_transit_idxs.append(mask_transit_idxs)

		##### update!
		self.mask_transit_idxs = np.array(all_quarter_mask_transit_idxs)

		if skip_quarter == 'n':
			try:
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
					try:
						detrend_model, fluxes_detrend, errors_detrend = cofiam_detrend(times=dtimes, fluxes=dfluxes, errors=derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)
					except:
						detrend_model, fluxes_detrend, errors_detrend = cofiam_detrend(times=np.array(dtimes, dtype=np.float64), fluxes=np.array(dfluxes, dtype=np.float64), errors=np.array(derrors, dtype=np.float64), telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)

					flags_detrend = dflags





				elif dmeth == 'polyAM':
					"""
					Polynomial detrending is the same basic algorithm as CoFiAM, but instead of using sinusoids, 
					it uses polynomials, and minimizes autocorrelation.
					"""
					max_degree = max_order(dtimes, self.duration_days)
					try:
						detrend_model, fluxes_detrend, errors_detrend = polyAM_detrend(times=dtimes, fluxes=dfluxes, errors=derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)
					except:
						detrend_model, fluxes_detrend, errors_detrend = polyAM_detrend(times=np.array(dtimes, type=np.float64), fluxes=np.array(dfluxes, type=np.float64), errors=np.array(derrors, type=np.float64), telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)

					flags_detrend = dflags





				elif dmeth == 'polyLOC':
					"""
					Like polyAM, except that it's a LOCAL fit that minimizes the BIC... AND it does one transit at a time, 
					so you need to handle this a bit differently

					"""

					detrend_model, fluxes_detrend, errors_detrend = [], [], []
					all_detrend_models, all_transit_times_detrend, all_transit_fluxes_detrend, all_transit_errors_detrend = [], [], [], []
					

					for ntau, tau in enumerate(self.taus):
						if tau < np.nanmin(dtimes) or tau > np.nanmax(dtimes):
							continue 

						print('tau '+str(ntau)+' of '+str(len(self.taus)))
						#### FOR EVERY TRANSIT!
						local_window_duration = 10*self.duration_days
						lwdm = 10 ### "local window duration multiple"
						if local_window_duration > 0.5 * self.period:
							#### you need a shorter duration!
							local_window_duration = 0.5 * self.period
							lwdm = local_window_duration / self.duration_days
						###local_window_duration_multiple
						
						max_degree = max_order(local_window_duration, self.duration_days)

						detrend_model, transit_times_detrend, transit_fluxes_detrend, transit_errors_detrend = polyLOC_detrend(times=np.array(dtimes, dtype=np.float64), fluxes=np.array(dfluxes, dtype=np.float64), errors=np.array(derrors, dtype=np.float64), local_window_duration_multiple=lwdm, telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)
						
						#### APPEND THIS LIST OF VALUES OT A BIGGER LIST
						all_detrend_models.append(detrend_model)
						all_transit_times_detrend.append(transit_times_detrend)
						all_transit_fluxes_detrend.append(transit_fluxes_detrend)
						all_transit_errors_detrend.append(transit_errors_detrend)

					### MAKE THE BIG LIST A BIG ARRAY
					try:
						all_detrend_models = np.concatenate(all_detrend_models)
						all_transit_times_detrend = np.concatenate(all_transit_times_detrend)
						all_transit_fluxes_detrend = np.concatenate(all_transit_fluxes_detrend)
						all_transit_errors_detrend = np.concatenate(all_transit_errors_detrend)
					except:
						print('could not concatenate (maybe an empty quarter')

					#### now we need to replace all non-transiting fluxes_detrend and errors_detrend with NaNs
					##### WHY? SO THAT IT IS THE SAME SIZE AS times, fluxes, and errors (for the way it saves in the file).
					for dtidx, dt in enumerate(dtimes):
						if dt in all_transit_times_detrend:
							detrend_model.append(all_detrend_models[dtidx])
							fluxes_detrend.append(all_transit_fluxes_detrend[dtidx])
							errors_detrend.append(all_transit_errors_detrend[dtidx])
						elif dt not in all_transit_times_detrend:
							detrend_model.append(np.nan)
							fluxes_detrend.append(np.nan)
							errors_detrend.append(np.nan)
						flags_detrend = dflags 




				elif dmeth == 'untrendy':
					print("UNTRENDY IMPLEMENTATION IS STILL BUGGY. BEWARE!")
					mask_transit_idxs = mask_transit_idxs.astype(int)
					print('mask_transit_idxs = ', mask_transit_idxs)
					detrend_model, fluxes_detrend, errors_detrend = untrendy_detrend(times=dtimes, fluxes=dfluxes, errors=derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs)
					#flags_detrend = np.linspace(0,0,len(fluxes_detrend))
					flags_detrend = dflags



				elif dmeth == 'george':
					print("GEORGE IMPLEMENTATION IS STILL BUGGY. BEWARE!")
					detrend_model, fluxes_detrend, errors_detrend = george_detrend(times=dtimes, fluxes=dfluxes, errors=derrors, GP_kernel=GP_kernel, metric=GP_metric, telescope=self.telescope, mask_idxs=mask_transit_idxs)
					#flags_detrend = np.linspace(0,0,len(fluxes_detrend))
					flags_detrend = dflags



				elif dmeth == 'medfilt':
					print("MEDIAN FILTERING HAS YET TO BE TESTED EXTENSIVELY. BEWARE!")
					detrend_model, fluxes_detrend, errors_detrend = medfilt_detrend(times=dtimes, fluxes=dfluxes, errors=derrors,  kernel_hours=(medfilt_kernel_transit_multiple*self.duration_hours), telescope=self.telescope, mask_idxs=mask_transit_idxs)
					flags_detrend = dflags


				elif dmeth == 'methmarg':
					local_window_duration = 10*self.duration_days
					lwdm = 10 ### "local window duration multiple"
					if local_window_duration > 0.5 * self.period:
						#### you need a shorter duration!
						local_window_duration = 0.5 * self.period
						lwdm = local_window_duration / self.duration_days

					
					print('RUNNING METHOD MARGINALIZATION! MAKE SURE YOU CHECK YOUR INPUTS!')
					detrend_model, fluxes_detrend, errors_detrend = methmarg_detrend(times=dtimes, fluxes=dfluxes, errors=derrors,  GP_kernel=GP_kernel, metric=GP_metric, kernel_hours=(medfilt_kernel_transit_multiple*self.duration_hours), local_window_duration_multiple=lwdm, telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)
					flags_detrend = dflags



			except:
				traceback.print_exc()
				print('DETRENDING FAILED FOR THIS QUARTER.')
				detrend_model, fluxes_detrend, errors_detrend = dfluxes, dfluxes, derrors 
				flags_detrend = np.linspace(2097152,2097152,len(fluxes_detrend))
				exceptions_raised = 'y'


		### update self -- just this quarter!
		if (skip_ntqs == 'n') and (exceptions_raised == 'n'):
			assert np.all(dfluxes != fluxes_detrend)
			assert np.all(derrors != errors_detrend)


		if len(self.quarters) > 1:
			master_detrend_model.append(np.array(detrend_model))
			master_detrend.append(np.array(fluxes_detrend))
			master_error_detrend.append(np.array(errors_detrend))
			master_flags_detrend.append(np.array(flags_detrend))
		elif len(self.quarters) == 1:
			master_detrend_model = np.array(master_detrend_model)
			master_detrend = np.array(fluxes_detrend)
			master_error_detrend = np.array(errors_detrend)
			master_flags_detrend = np.array(flags_detrend)


	### this is the first initialization of the detrended fluxes.
	self.detrend_model = master_detrend_model
	self.fluxes_detrend = master_detrend
	self.errors_detrend = master_error_detrend
	self.flags_detrend = master_flags_detrend 

	### before you overwrite the flags, compare them.

	if len(self.quarters) == 1:
		final_flags = np.nanmax((self.flags, self.flags_detrend), axis=0)

	else:
		final_flags = []
		for qidx in np.arange(0,len(self.quarters),1):
			qfinal_flags = np.nanmax((self.flags[qidx], self.flags_detrend[qidx]), axis=0)
			final_flags.append(np.array(qfinal_flags))
		self.flags_detrend = final_flags

	if save_lc == 'y':
		lcfile = open(self.savepath+'/'+self.target+'_'+self.telescope+'_lightcurve.tsv', mode='w')
		if self.telescope.lower() == 'kepler' or self.telescope.lower() == 'k2':
			lcfile.write('BKJD\tfluxes\terrors\tdetrend_model\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
		elif self.telescope.lower() == 'tess':
			lcfile.write('BTJD\tfluxes\terrors\tdetrend_model\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
		elif self.telescope.lower() == 'user':
			lcfile.write('BJD\tfluxes\terrors\tdetrend_model\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
		
		### overwrite the existing file!
		self.detrend_model = np.array(self.detrend_model, dtype=object)
		self.fluxes_detrend = np.array(self.fluxes_detrend, dtype=object)
		self.errors_detrend = np.array(self.errors_detrend, dtype=object)
		self.flags_detrend = np.array(self.flags_detrend, dtype=object)
		lc_times, lc_fluxes, lc_errors, lc_detrend_model, lc_fluxes_detrend, lc_errors_detrend, lc_flags = self.times, self.fluxes, self.errors, self.detrend_model, self.fluxes_detrend, self.errors_detrend, self.flags_detrend

		if len(self.quarters) > 1:
			for qidx in np.arange(0,len(self.quarters),1):
				qtq = self.quarters[qidx]
				qtimes, qfluxes, qerrors, qdetrend_model, qfluxes_detrend, qerrors_detrend, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_detrend_model[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx]

				for qt, qf, qe, qdm, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qdetrend_model, qfluxes_detrend, qerrors_detrend, qflags):
					lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qdm)+'\t'+str(qfd)+'\t'+str(qed)+'\t'+str(qfl)+'\t'+str(qtq)+'\n')

				print('quarter written to file.')

		elif len(self.quarters) == 1:
			qtimes, qfluxes, qerrors, qdetrend_model, qfluxes_detrend, qerrors_detrend, qflags = lc_times, lc_fluxes, lc_errors, lc_detrend_model, lc_fluxes_detrend, lc_errors_detrend, lc_flags
			
			if qtimes.shape[0] == 1:
				qtimes = qtimes[0]
			if qfluxes.shape[0] == 1:
				qfluxes = qfluxes[0]
			if qerrors.shape == 1:
				qerrors = qerrors[0]

			try:
				if qdetrend_model.shape == 1:
					qdetrend_model = qdetrend_model[0]
			except:
				if type(qdetrend_model) == list:
					qdetrend_model = qdetrend_model[0]

			try:
				if qfluxes_detrend.shape == 1:
					qfluxes_detrend = qfluxes_detrend[0]
			except:
				if type(qfluxes_detrend) == list:
					qfluxes_detrend = qfluxes_detrend[0]

			try:
				if qerrors_detrend.shape == 1:
					qerrors_detrend = qerrors_detrend[0]
			except:
				if type(qerrors_detrend) == list:
					qerrors_detrend = qerrors_detrend[0]

			if len(qflags) != len(qerrors_detrend):
				try:
					if qflags.shape == 1:
						qflags = qflags[0]
				except:
					if type(qflags) == list:
						qflags = qflags[0]

			for qt, qf, qe, qdm, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qdetrend_model, qfluxes_detrend, qerrors_detrend, qflags):
				lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qfd)+'\t'+str(qed)+'\t'+str(qfl)+'\t'+str(self.quarters[0])+'\n')

		lcfile.close()

	if self.telescope.lower() != 'user':
		self.get_neighbors(save_to_file='y')

	### finally, run get_neighbors() so that the neighbors will be appended to the end of the file.





def gen_batman(self, folded='n'):
	from mp_batman import run_batman 

	"""
	##### this method will generate a batman light curve for the light curve object.
	#run_batman(all_times, RpRstar, Rstar, bplan, Pplan, tau0, q1, q2, long_peri=0, ecc=0, Mstar=None, 
		Mplan=None, rhostar=None, rhoplan=None, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', 
		ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n', **kwargs):
	"""
	if np.isfinite(self.longp):
		bat_longp = self.longp
	else:
		bat_longp = 0

	if np.isfinite(self.eccen):
		bat_eccen = self.eccen
	else:
		bat_eccen = 0

	bat_sma = self.sma_AU * au.value 

	#if folded == 'n':
	bat_times, bat_fluxes = run_batman(all_times=self.times, RpRstar=self.rprstar, Rstar=self.rstar_meters, bplan=self.impact, Pplan=self.period, tau0=self.tau0, q1=self.q1, q2=self.q2, planet_sma=bat_sma)
	self.bat_times = bat_times
	self.bat_fluxes = bat_fluxes 	

	if folded == 'y':
		#sorted_times_idxs = np.argsort(times)
		#folded_bat_times, folded_bat_fluxes = run_batman(all_times=times, RpRstar=self.rprstar, Rstar=self.rstar_meters, bplan=self.impact, Pplan=self.period, tau0=0, q1=self.q1, q2=self.q2, planet_sma=bat_sma)
		folded_bat_times, folded_bat_fluxes, folded_bat_errors = lc_fold(times=self.bat_times, fluxes=self.bat_fluxes, errors=self.errors_detrend, tau0=self.tau0, period=self.period)

		self.folded_bat_times = folded_bat_times
		self.folded_bat_fluxes = folded_bat_fluxes



