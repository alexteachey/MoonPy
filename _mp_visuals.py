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
from matplotlib.offsetbox import AnchoredText

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
from mp_animate import * 
try:
	import pandoramoon as pandora 
	from pandoramoon.helpers import ld_convert, ld_invert 
except:
	print('Pandora functions did not load. Maybe not installed.')


plt.rcParams["font.family"] = 'serif'

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/_mp_visuals.py')]


def plot_lc(self, facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', folded='n', include_flagged='n', undetrended='y', detrended='y', show_errors='n', show_stats='y', show_neighbors='y', mask_multiple=None, period=None, show_model='y', show_batman='y', show_pandora='y', show_model_residuals='y', time_format='native', pltshow='y', phase_offset=0.0, binned='n'):
	print('calling _mp_visuals.py/plot_lc().')
	#if ('detrend_model' not in dir(self)) or (np.any(np.isfinite(np.concatenate(self.detrend_model))) == False):
	#	detrended = 'n'
	#	show_model = 'n'
	#### CHETAN'S IMPROVEMENT vvv 
	if len(self.quarters) > 1:
	    if ('detrend_model' not in dir(self)) or (np.any(np.isfinite(np.concatenate(self.detrend_model))) == False):
	        detrended = 'n'
	        show_model = 'n'
	elif len(self.quarters) == 1:
	    if ('detrend_model' not in dir(self)) or (np.any(np.isfinite(np.array(self.detrend_model, dtype=np.float64))) == False):
	        detrended = 'n'
	        show_model = 'n'


	if self.telescope.lower() == 'k2' and include_flagged == 'n':
		print(' ')
		print('BE ADVISED: K2 photometry may have lots of excluded data points due to flagging.')
		print("call plot_lc(include_flagged='y') if the light curve is sparse.")
		print(' ')

	if period == None:
		try:
			period = self.period 
		except:
			period = np.nan

	### THIS FUNCTION PLOTS THE LIGHT CURVE OBJECT.
	
	#if mask_multiple == None:
	#	mask_multiple = self.mask_multiple

	if detrended == 'n':
		undetrended = 'y'

	if (undetrended == 'y') and (detrended == 'y'):
		nplots = 2
		if folded=='n':
			fig, ax = plt.subplots(2, sharex=True, figsize=(6,8))
		elif folded=='y':
			fig, ax = plt.subplots(2, figsize=(6,8))

	else:
		#### will be just one or the other
		nplots = 1
		fig, ax = plt.subplots()

	try:
		plot_times, plot_fluxes, plot_errors, plot_fluxes_detrend, plot_errors_detrend, plot_flags, plot_quarters = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags, self.quarters

	except:
		print("WARNING: light curve has not been detrended yet!")
		detrended = 'n'
		nplots = 1
		undetrended = 'y'
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
		stitched_times, stitched_fluxes, stitched_errors, stitched_fluxes_detrend, stitched_errors_detrend, stitched_flags, stitched_quarters = np.hstack((plot_times)), np.hstack((plot_fluxes)), np.hstack((plot_errors)), np.hstack((plot_fluxes_detrend)), np.hstack((plot_errors_detrend)), np.hstack((plot_flags)), np.hstack((plot_quarters))
	else:
		stitched_times, stitched_fluxes, stitched_errors, stitched_flags, stitched_quarters = np.hstack((plot_times)), np.hstack((plot_fluxes)), np.hstack((plot_errors)), np.hstack((plot_flags)), np.hstack((plot_quarters))			

	if include_flagged=='n':
		### remove all data points with qflag != 0
		badflag_idxs = np.where(stitched_flags != 0)[0]
		if detrended == 'y':
			stitched_times, stitched_fluxes, stitched_errors, stitched_fluxes_detrend, stitched_errors_detrend, stitched_flags = np.delete(stitched_times, badflag_idxs), np.delete(stitched_fluxes, badflag_idxs), np.delete(stitched_errors, badflag_idxs), np.delete(stitched_fluxes_detrend, badflag_idxs), np.delete(stitched_errors_detrend, badflag_idxs), np.delete(stitched_flags, badflag_idxs)
		else:
			stitched_times, stitched_fluxes, stitched_errors, stitched_flags = np.delete(stitched_times, badflag_idxs), np.delete(stitched_fluxes, badflag_idxs), np.delete(stitched_errors, badflag_idxs), np.delete(stitched_flags, badflag_idxs)

		assert np.all(stitched_flags == 0)


	#### IDENTIFY TARGET TIMES
	try:
		target_taus = self.taus 
		target_dur = self.duration_days
	except:
		target_taus = np.array([np.nan])
		target_dur = np.array([np.nan]) 

	target_transit_idxs = []
	for tt in target_taus:
		ttidxs = np.where((stitched_times >= (tt - (self.mask_multiple/2)*target_dur)) & (stitched_times <= (tt + (self.mask_multiple/2)*target_dur)))[0]
		target_transit_idxs.append(ttidxs)
	target_transit_idxs = np.hstack(target_transit_idxs)


	### this will highlight all the other transits for the neighbors (if any)
	try:
		neighbors = self.neighbor_dict.keys()
	except:
		print('NEIGHBOR DICTIONARY UNAVAILABLE.')
		self.neighbor_dict = {}
		self.neighbors = []
		neighbors = np.array([])

	if time_format == 'native':
		plot_stitched_times = stitched_times
	elif time_format == 'bjd':
		if self.telescope == 'kepler':
			plot_stitched_times = stitched_times + 2454833
		elif self.telescope == 'tess':
			plot_stitched_times = stitched_times + 2457000

	if folded == 'n':

		if nplots == 2:
			ax[0].scatter(plot_stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			ax[1].scatter(plot_stitched_times, stitched_fluxes_detrend, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			ax[0].set_xlim(np.nanmin(plot_stitched_times), np.nanmax(plot_stitched_times))
			ax[1].set_xlim(np.nanmin(plot_stitched_times), np.nanmax(plot_stitched_times))			

			if show_errors == 'y':
				ax[0].errorbar(plot_stitched_times, stitched_fluxes, yerr=stitched_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
				ax[1].errorbar(plot_stitched_times, stitched_fluxes_detrend, yerr=stitched_errors_detrend, ecolor='k', zorder=0, alpha=0.5, fmt='none')

			if show_model == 'y':
				try:
					#ax[0].plot(plot_stitched_times, np.concatenate(self.detrend_model), color='BlueViolet', linewidth=2, alpha=0.7)
	                ##### CHETAN'S IMRPOVEMENT vvv
					if (quarters=='all') and (len(self.quarters) > 1):
						ax[0].plot(plot_stitched_times, np.concatenate(self.detrend_model), color='BlueViolet', linewidth=2, alpha=0.7)
					elif (quarters=='all') and (len(self.quarters) == 1):
						ax[0].plot(plot_stitched_times, np.array(self.detrend_model, dtype=np.float64), color='BlueViolet', linewidth=2, alpha=0.7)
					elif (quarters!='all') and (len(quarters) > 1):
						ax[0].plot(plot_stitched_times, np.concatenate(self.detrend_model), color='BlueViolet', linewidth=2, alpha=0.7)
					elif (quarters!='all') and (len(quarters) == 1):
						ax[0].plot(plot_stitched_times, np.array(self.detrend_model, dtype=np.float64), color='BlueViolet', linewidth=2, alpha=0.7)

				except:
					print('self.detrend_model not stored. (FIX THIS BUG).')

		elif nplots == 1: ### detrended or undetrended, but not both
			if detrended == 'y':
				ax.scatter(plot_stitched_times, stitched_fluxes_detrend, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
				ax.set_xlim(np.nanmin(plot_stitched_times), np.nanmax(plot_stitched_times))							
				if show_errors == 'y':
					ax.errorbar(plot_stitched_times, stitched_fluxes_detrend, yerr=stitched_errors_detrend, ecolor='k', zorder=0, alpha=0.5, fmt='none')

			else:
				ax.scatter(plot_stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
				if show_errors == 'y':
					ax.errorbar(plot_stitched_times, stitched_fluxes, yerr=stitched_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
				if show_model == 'y':
					try:
						#ax.plot(plot_stitched_times, np.concatenate(self.detrend_model), color='BlueViolet', linewidth=2, alpha=0.7)
						#### CHETAN'S IMPROVEMENT vvvv
						if (quarters=='all') and (len(self.quarters) > 1):
							ax.plot(plot_stitched_times, np.concatenate(self.detrend_model), color='BlueViolet', linewidth=2, alpha=0.7)
						elif (quarters=='all') and (len(self.quarters) == 1):
							ax.plot(plot_stitched_times, np.array(self.detrend_model, dtype=np.float64), color='BlueViolet', linewidth=2, alpha=0.7)
						elif (quarters!='all') and (len(quarters) > 1):
							ax.plot(plot_stitched_times, np.concatenate(self.detrend_model), color='BlueViolet', linewidth=2, alpha=0.7)
						elif (quarters!='all') and (len(quarters) == 1):
							ax.plot(plot_stitched_times, np.array(self.detrend_model, dtype=np.float64), color='BlueViolet', linewidth=2, alpha=0.7)

					except:
						print('detrend_model not available.')


		for neighbor in neighbors:
			try:
				neighbor_taus = self.neighbor_dict[neighbor].taus 
			except:
				neighbor_taus = np.array([])

			try:
				neighbor_dur = self.neighbor_dict[neighbor].duration_days 
			except:
				neighbor_dur = np.array([])

			try:
				neighbor_transit_idxs = []
				for nt in neighbor_taus:
					ntidxs = np.where((stitched_times >= (nt - (self.mask_multiple/2)*neighbor_dur)) & (stitched_times <= (nt + (self.mask_multiple/2)*neighbor_dur)))[0]
					neighbor_transit_idxs.append(ntidxs)
				neighbor_transit_idxs = np.hstack(neighbor_transit_idxs)
				
				if (nplots == 2) and (show_neighbors == 'y'):
					ax[0].scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
					ax[1].scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes_detrend[neighbor_transit_idxs], s=10, marker='x', label=neighbor)

				elif (nplots == 1) and (show_neighbors == 'y'):
					if detrended == 'y':
						ax.scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes_detrend[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
					else:
						ax.scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
			except:
				traceback.print_exc()

		### PLOT THE TARGET TRANSITS TOO!
		if (nplots == 2) and (show_neighbors == 'y'):
			ax[0].scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')
			ax[1].scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes_detrend[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')			

		elif (nplots == 1) and (show_neighbors == 'y'):
			if detrended == 'y':
				ax.scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes_detrend[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')	
			else:
				ax.scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')



		if (show_batman == 'y') and (detrended == 'y'):
			try:
				self.gen_batman(folded='n')
				if nplots == 2:
					ax[1].plot(self.bat_times[qstokeep_idxs], self.bat_fluxes[qstokeep_idxs], c='BlueViolet', linewidth=2, zorder=5, alpha=0.7, label='planet model')	
				elif nplots == 1:
					ax.plot(self.bat_times[qstokeep_idxs], self.bat_fluxes[qstokeep_idxs], c='BlueViolet', linewidth=2, zorder=5, alpha=0.7, label='planet model')	

					
			except:
				print("COULD NOT GENERATE A BATMAN MODEL FOR THIS PLANET.")


		#### NEW AUGUST 2022 -- show a pandora model! 
		if (show_pandora == 'y') and (detrended == 'y'):

			try:
				self.get_Pandora_posteriors(model='M')
			except:
				print('unable to get Pandora posteriors for model M.')

			try:
				self.get_Pandora_posteriors(model='P')
			except:
				print('unable to get Pandora posteriors for model P.')

			models_available = []

			try:
				print('Planet posteriors available: ')
				print(self.Pandora_planet_PEWdict.keys())
				planet_only_avail = True 
				models_available.append('P')
			except:
				print('Planet posteriors not available.')
				planet_only_avail = False 

			try:
				print('Moon posteriors available: ')
				print(self.Pandora_moon_PEWdict.keys())
				planet_plus_moon_avail = True
				models_available.append('M')
			except:
				print('Moon posteriors not available.')
				planet_plus_moon_avail = False 

			if len(models_available) > 0:


				for model in models_available:
					
					if model.lower() == 'p':
						model_PEWdict = self.Pandora_planet_PEWdict 
					elif model.lower() == 'm':
						model_PEWdict = self.Pandora_moon_PEWdict 

				ndraws = 100
				nsamples = len(model_PEWdict['q1'])
				random_idxs = np.random.choice(np.arange(0,nsamples,1), size=ndraws)


				#### compute median values
				for key in model_PEWdict.keys():
					print('median '+str(key)+': '+str(np.nanmedian(model_PEWdict[key])))



				##### generate a model from the NASA Exoplanet Archive
				#### other fixed parameters

				#### conversion 

				"""
				medq1, medq2 = np.nanmedian(model_PEWdict['q1']), np.nanmedian(model_PEWdict['q2'])
				medu1, medu2 = ld_convert(medq1, medq2) 

				#### NOW GENERATE THE MODEL! 	
				NEAparams = pandora.model_params()
				
				NEAparams.R_star = self.st_rad * R_sun.value 
				NEAparams.per_bary = self.pl_orbper 
				NEAparams.a_bary = (self.pl_orbsmax * au.value) / (self.st_rad * R_sun.value) # a/Rstar
				NEAparams.r_planet = (self.pl_rade * R_earth.value) / (self.st_rad * R_sun.value) #Rp/Rstar
				NEAparams.b_bary = self.pl_imppar
				NEAparams.t0_bary_offset = np.nanmedian(model_PEWdict['t0_bary_offset']) ### this is not in NEA 
				NEAparams.M_planet = self.pl_bmasse * M_earth.value # kg
				NEAparams.ecc_bary = self.pl_orbeccen
				NEAparams.w_bary = 83
				#### moon parameters
				NEAparams.r_moon = 1e-8 #### PARAM #7 -- satellite radius divided by stellar radius
				NEAparams.per_moon = 30 #### PARAM #8 -- need to define above
				NEAparams.tau_moon = 0 #### PARAM #9-- must be between zero and one -- I think this is the phase...
				NEAparams.Omega_moon = 0 #### PARAM # 10 -- longitude of the ascending node??? between 0 and 180 degrees 
				NEAparams.i_moon = 0 #### PARAM #11 -- between 0 and 180 degrees
				NEAparams.M_moon = 1e-8 
				NEAparams.u1 = float(medu1) #### PARAM #13 need to define above!
				NEAparams.u2 = float(medu2) #### PARAM #14 need to define above!
				NEAparams.t0_bary = self.tau0 
				NEAparams.epochs = len(self.taus) #### needs to be defined above!
				NEAparams.epoch_duration = 3 #### need to be defined above
				NEAparams.cadences_per_day = 48 
				#params.epoch_distance = Pplan 
				NEAparams.epoch_distance = self.pl_orbper
				NEAparams.supersampling_factor = 1
				NEAparams.occult_small_threshold = 0.1 ### between 0 and 1 -- what is this?
				NEAparams.hill_sphere_threshold = 1.2 #### what does this mean?

				NEApdtime = pandora.time(NEAparams).grid()
				NEApdmodel = pandora.moon_model(NEAparams)

				total_flux, planet_flux, moon_flux = NEApdmodel.light_curve(NEApdtime)


				if nplots == 2:
					ax[1].plot(NEApdtime, total_flux, c='green', linewidth=2, zorder=100, alpha=0.7, label='NEA')	
				elif nplots == 1:
					ax.plot(NEApdtime, total_flux, c='green', linewidth=2, zorder=100, alpha=0.7, label='NEA')						

				"""


				#### NOW STEP THROUGH THE POSTERIORS
				for nidx, random_idx in enumerate(random_idxs):
					q1 = model_PEWdict['q1'][random_idx]
					q2 = model_PEWdict['q2'][random_idx]
					per_bary = model_PEWdict['per_bary'][random_idx]
					a_bary = model_PEWdict['a_bary'][random_idx]
					r_planet = model_PEWdict['r_planet'][random_idx]
					b_bary = model_PEWdict['b_bary'][random_idx]
					w_bary = model_PEWdict['w_bary'][random_idx]
					ecc_bary = model_PEWdict['ecc_bary'][random_idx]
					t0_bary_offset = model_PEWdict['t0_bary_offset'][random_idx]

					if model.lower() == 'p':
						#### set the moon values to standard values for no moon present
						r_moon = 1e-8
						per_moon = 30 
						tau_moon = 0
						Omega_moon = 0
						i_moon = 0
						ecc_moon = 0
						w_moon = 0
						M_moon = 1e-8 

					elif model.lower() == 'm':
						#### uses the posteriors! 
						r_moon = model_PEWdict['r_moon'][random_idx]
						per_moon = model_PEWdict['per_moon'][random_idx]
						tau_moon = model_PEWdict['tau_moon'][random_idx]
						Omega_moon = model_PEWdict['Omega_moon'][random_idx]
						i_moon = model_PEWdict['i_moon'][random_idx]
						ecc_moon = model_PEWdict['ecc_moon'][random_idx]
						w_moon = model_PEWdict['w_moon'][random_idx]
						M_moon = model_PEWdict['M_moon'][random_idx]


					#### other fixed parameters
					R_star = self.st_rad * R_sun.value 
					t0_bary = self.tau0 
					M_planet = self.pl_bmasse * M_earth.value

					#### conversion 
					u1, u2 = ld_convert(q1, q2) 

					#### NOW GENERATE THE MODEL! 
					params = pandora.model_params()
					
					params.R_star = float(R_star) #### FIT PARAM #0 
					params.per_bary = float(per_bary) #### PARAM #1 Pplan [days]
					params.a_bary = float(a_bary)  #### PARAM #2 
					params.r_planet = float(r_planet) #### PARAM #3 
					params.b_bary = float(b_bary) #### PARAM # 4
					params.t0_bary_offset = float(t0_bary_offset) #### PARAM #5 what is this? 
					params.M_planet = float(M_planet) #### PARAM #6 [kg]
					params.ecc_bary = float(ecc_bary)
					params.w_bary = float(w_bary)
					#### moon parameters
					params.r_moon = float(r_moon) #### PARAM #7 -- satellite radius divided by stellar radius
					params.per_moon = float(per_moon) #### PARAM #8 -- need to define above
					params.tau_moon = float(tau_moon) #### PARAM #9-- must be between zero and one -- I think this is the phase...
					params.Omega_moon = float(Omega_moon) #### PARAM # 10 -- longitude of the ascending node??? between 0 and 180 degrees 
					params.i_moon = float(i_moon) #### PARAM #11 -- between 0 and 180 degrees
					params.M_moon = float(M_moon) 
					params.u1 = float(u1) #### PARAM #13 need to define above!
					params.u2 = float(u2) #### PARAM #14 need to define above!


					#### FIXED PARAMETERS -- BUT I"m NOT SURE WHY ... 
					params.t0_bary = float(t0_bary) #### FIXED?

					#### other inputs
					params.epochs = len(self.taus) #### needs to be defined above!
					params.epoch_duration = 3 #### need to be defined above
					params.cadences_per_day = 48 
					#params.epoch_distance = Pplan 
					params.epoch_distance = per_bary 
					params.supersampling_factor = 1
					params.occult_small_threshold = 0.1 ### between 0 and 1 -- what is this?
					params.hill_sphere_threshold = 1.2 #### what does this mean?

					pdtime = pandora.time(params).grid()
					pdmodel = pandora.moon_model(params)

					#total_flux, planet_flux, moon_flux = pdmodel.light_curve(pdtime)
					#total_flux, planet_flux, moon_flux = pdmodel.light_curve(all_times)
					total_flux, planet_flux, moon_flux = pdmodel.light_curve(np.concatenate(self.times))

					#### plot them!
					if model.lower() == 'p':
						model_label = 'Pandora (planet only)'
						model_color = 'DarkOrange'

					elif model.lower() == 'm':
						model_label = 'Pandora (planet+moon)'
						model_color = 'BlueViolet'

					if ndraws > 1:
						linewidth=1
						alpha=0.3
					elif ndraws == 1:
						linewidth=2
						alpha=0.7

					if nidx == 0:
						if nplots == 2:
							ax[1].plot(np.concatenate(self.times), total_flux, c=model_color, linewidth=linewidth, zorder=5, alpha=alpha, label=model_label)	
						elif nplots == 1:
							ax.plot(np.concatenate(self.times), total_flux, c=model_color, linewidth=linewidth, zorder=5, alpha=alpha, label=model_label)	

					else:
						#### don't label
						if nplots == 2:
							ax[1].plot(np.concatenate(self.times), total_flux, c=model_color, linewidth=linewidth, zorder=5, alpha=alpha)	
						elif nplots == 1:
							ax.plot(np.concatenate(self.times), total_flux, c=model_color, linewidth=linewidth, zorder=5, alpha=alpha)	


					"""
					if print_params == 'y':

						print(" ")
						print("Rp/Rstar = ", RpRstar)
						print("transit depth [ppm] = ", RpRstar**2 * 1e6)
						print("stellar density [kg / m^3] = ", rhostar)
						print("impact = ", bplan)
						print("Period [days] = ", Pplan)
						print("tau_0 [day] = ", tau0)
						print("q1,q2 = ", q1, q2)
						if (model == 'M') or (model == "Z"):
							print("planet density [kg / m^3] = ", rhoplan)
							print("sat_sma = [Rp] ", sat_sma)
							print("sat_phase = ", sat_phase)
							print("sat_inc = ", sat_inc)
							print("sat_omega = ", sat_omega)
							print("Msat / Mp = ", MsatMp)
							print("Rsat / Rp = ", RsatRp)
						print(" ")
					"""


					#return output_times, output_fluxes 
					#return total_flux, planet_flux, moon_flux





	elif folded == 'y':
		nplots = 1 #### should only show the detrend -- folding on undetrended is nonsense.
		plt.close()
		fig, ax = plt.subplots()
		detrended = 'y' #### it doesn't make any sense to phase-fold an undetrended light curve
		try:
			self.fold(detrended=detrended, phase_offset=phase_offset, period=period)
		except:
			self.get_properties(locate_neighbor='n')
			self.fold(phase_offset=phase_offset, period=period)


		#### PLOT THE NEIGHBORS 
		for neighbor in neighbors:
			neighbor_taus = self.neighbor_dict[neighbor].taus 
			neighbor_dur = self.neighbor_dict[neighbor].duration_days 

			neighbor_transit_idxs = []
			for nt in neighbor_taus:
				ntidxs = np.where((plot_stitched_times >= (nt - (self.mask_multiple/2)*neighbor_dur)) & (plot_stitched_times <= (nt + (self.mask_multiple/2)*neighbor_dur)))[0]
				neighbor_transit_idxs.append(ntidxs)
			neighbor_transit_idxs = np.hstack(neighbor_transit_idxs)
			
			if (nplots == 2) and (show_neighbors == 'y'):
				ax[0].scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
				#### DOESN'T MAKE SENSE TO PLOT NEIGHBORS IN THE PHASE FOLD
				#ax[1].scatter(self.fold_times[neighbor_transit_idxs], self.fold_fluxes[neighbor_transit_idxs], s=10, marker='x', label=neighbor)

			#### DOESN'T MAKE SENSE TO SHOW NEIGHBORS IN THE PHASE FOLD!
			#elif (nplots == 1) and (show_neighbors == 'y'):
				#if detrended == 'y':
				#	ax.scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes_detrend[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
				#else:
				#	ax.scatter(self.fold_times[neighbor_transit_idxs], self.fold_fluxes[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
		

		### PLOT THE TARGET TRANSITS TOO!
		if (nplots == 2) and (show_neighbors == 'y'):
			ax[0].scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')
			#ax[1].scatter(self.fold_times[target_transit_idxs], self.fold_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')			

		#elif (nplots == 1) and (show_neighbors == 'y'):
			#if detrended == 'y':
			#	ax.scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes_detrend[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')	
			#else:
			#	ax.scatter(self.fold_times[target_transit_idxs], self.fold_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')


		if binned == 'n':	
			if nplots == 2:
				ax[0].scatter(plot_stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
				ax[1].scatter(self.fold_times, self.fold_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
				if show_errors == 'y':
					ax[0].errorbar(plot_stitched_times, stitched_fluxes, yerr=stitched_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
					ax[1].errorbar(self.fold_times, self.fold_fluxes, yerr=self.fold_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
					ax[0].set_xlim(np.nanmin(plot_stitched_times), np.nanmax(plot_stitched_times))
					ax[1].set_xlim(np.nanmin(self.fold_times), np.nanmax(self.fold_times))		

			elif nplots == 1:
				ax.scatter(self.fold_times, self.fold_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
				if show_errors == 'y':
					ax.errorbar(self.fold_times, self.fold_fluxes, yerr=self.fold_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
					ax.set_xlim(np.nanmin(self.fold_times), np.nanmax(self.fold_times))		

		elif binned == 'y':
			detrended = 'y'
			undetrended = 'n'
			nplots = 1

			fold_bin_step = 0.0005
			fold_bins = np.arange(np.nanmin(self.fold_times), np.nanmax(self.fold_times), fold_bin_step)
			fold_bin_fluxes = []
			fold_bin_errors = []
			for fb in fold_bins:
				fb_idxs = np.where((self.fold_times >= fb- fold_bin_step/2) & (self.fold_times < fb + fold_bin_step/2))[0]
				fold_bin_fluxes.append(np.nanmedian(self.fold_fluxes[fb_idxs]))
				fold_bin_errors.append(np.nanstd(self.fold_fluxes[fb_idxs])/np.sqrt(len(fb_idxs)))

			fold_bin_fluxes, fold_bin_errors = np.array(fold_bin_fluxes), np.array(fold_bin_errors)

			ax.scatter(self.fold_times, self.fold_fluxes, facecolors='k', s=5, zorder=0, alpha=0.2)
			ax.scatter(fold_bins, fold_bin_fluxes, facecolor=facecolor, alpha=0.7, s=15, zorder=1)
			ax.set_xlim(np.nanmin(self.fold_times), np.nanmax(self.fold_times))

		if (show_batman == 'y') and (detrended == 'y'):
			try:
				self.gen_batman(folded='y')
				if nplots == 2:
					ax[1].plot(self.folded_bat_times, self.folded_bat_fluxes, c='BlueViolet', linewidth=2, zorder=5, alpha=0.7, label='planet model')	
					ax[0].set_xlim(np.nanmin(self.folded_bat_times), np.nanmax(self.folded_bat_times))
					ax[1].set_xlim(np.nanmin(self.folded_bat_times), np.nanmax(self.folded_bat_times))					
				elif nplots == 1:
					ax.plot(self.folded_bat_times, self.folded_bat_fluxes, c='BlueViolet', linewidth=2, zorder=5, alpha=0.7, label='planet model')	
					ax.set_xlim(np.nanmin(self.folded_bat_times), np.nanmax(self.folded_bat_times))					
					
			except:
				print("COULD NOT GENERATE A BATMAN MODEL FOR THIS PLANET.")





	if (self.telescope.lower() == 'kepler') or (self.telescope.lower() == 'k2'):
		if folded=='y':
			if nplots == 2:
				ax[1].set_xlabel('Phase')
			elif nplots == 1:
				ax.set_xlabel('Phase')
		else:
			if nplots == 2:
				ax[1].set_xlabel('BKJD')
			elif nplots == 1:
				ax.set_xlabel('BKJD')

	elif (self.telescope.lower() == 'tess'):
		if folded=='y':
			if nplots == 2:
				ax[1].set_xlabel('Phase')
			elif nplots == 1:
				ax.set_xlabel('Phase')
		else:
			if nplots == 2:
				ax[1].set_xlabel('BTJD')
			elif nplots == 1:
				ax.set_xlabel('BTJD')

	if nplots == 2:		
		ax[0].set_ylabel('Flux')
		ax[1].set_ylabel('Normalized Flux')
	elif nplots == 1:
		if detrended == 'y':
			ax.set_ylabel('Normalized Flux')
		elif detrended == 'n':
			ax.set_ylabel('Flux')


	try:
		if nplots == 2:
			ax[0].set_title(str(self.target))
		elif nplots == 1:
			ax.set_title(str(self.target))

	except:
		pass


	try:
		batman_transit_depth = (1 - np.nanmin(self.bat_fluxes))*1e6 #### ppm
	except:
		batman_transit_depth = np.nan


	##### ANALYZE WHETHER THE TARGET RESIDUALS ARE OFF
	try:
		full_LC_residuals = stitched_fluxes_detrend-self.bat_fluxes
		print('full_LC_residuals = ', full_LC_residuals)
		full_LC_median = np.nanmedian(full_LC_residuals)
		print('full_LC_median = ', full_LC_median)
		full_LC_std = np.nanstd(full_LC_residuals)
		print('full_LC_std = ', full_LC_std)
		full_LC_std_ppm = full_LC_std*1e6 
		print('full_LC_std_ppm = ', full_LC_std_ppm)

	except:
		full_LC_residuals = np.nan
		full_LC_median = np.nan 
		full_LC_std = np.nan 
		full_LC_std_ppm = np.nan 

	try:
		target_median = np.nanmedian(full_LC_residuals[target_transit_idxs])
		ntarget_points_outside_2sig = 0
		for ttp in full_LC_residuals[target_transit_idxs]:
			if (ttp > full_LC_median + 2*full_LC_std) or (ttp < full_LC_median - 2*full_LC_std):
				ntarget_points_outside_2sig += 1

		print('full_LC_median = ', full_LC_median)
		print('full_LC_std [ppm] = ', full_LC_std*1e6)
		fraction_target_points_outside_2sig = ntarget_points_outside_2sig / len(target_transit_idxs)
		print('fraction of target in-transit residuals outside 2sig: ', fraction_target_points_outside_2sig)
		if fraction_target_points_outside_2sig > 0.05:
			print("POSSIBLE BAD DETREND.")

	except:
		pass






	if (show_stats == 'y') and 'fluxes_detrend' in dir(self):
		textstr = '\n'.join((
		    r'depth $=%.2f$ ppm' % (batman_transit_depth, ),
		    r'scatter $=%.2f$ ppm' % (full_LC_std_ppm, ),))

		# these are matplotlib.patch.Patch properties
		props = dict(boxstyle='square', facecolor='white', alpha=0.5)

		# place a text box in upper left in axes coords
		#ax.text(0.75, 0.98, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
		#ax.text(0.65, 0.98, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
		anchored_text = AnchoredText(textstr, loc='lower left')
		if nplots == 2:
			ax[0].add_artist(anchored_text)
			ax[1].legend(loc='lower left')				
		elif nplots == 1:
			ax.add_artist(anchored_text)
			ax.legend(loc='upper right')

	if pltshow == 'y':	
		plt.show()
	else:
		pass


	if (show_model_residuals == 'y') and (show_batman == 'y') and (folded == 'n'):
		try:
			##### plot the light curve with the model removed
			plt.scatter(plot_stitched_times, stitched_fluxes_detrend-self.bat_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			plt.scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes_detrend[target_transit_idxs]-self.bat_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')		
			plt.xlabel('BKJD')
			plt.ylabel('fluxes - model')
			plt.show()

		except:
			print('BATMAN model not available.')




def animate_FFI(self, sector=None, ffi_idx=None, save_animation=False, clobber='n', normalize_flux=True, return_arrays=True):
	from mp_lcfind import TESS_direct_FFI_download 
	#def TESS_direct_FFI_download(tic, RA_hh=None, RA_deg=None, Dec_deg=None, npix_per_side=10, clobber='n'):
	TIC_filedir, FFI_files = TESS_direct_FFI_download(tic=self.target, RA_hh=self.RA, Dec_deg=self.Dec, clobber=clobber)

	for nFFI,FFI in enumerate(FFI_files):
		print(nFFI, FFI)

	if type(sector) == type(None) and type(ffi_idx) == type(None):
		#### as the user to specify which sector index they want to animate
		for nFFI,FFI in enumerate(FFI_files):
			print(nFFI, FFI)

		print(' ')
		ffi_idx = int(input('Enter the index number of the sector you want to animate: '))

	elif type(sector) != type(None):
		#### need to figure out the ffi_idx of this
		sector_string = str(sector)
		while len(sector_string) < 4:
			sector_string = '0'+sector_string
		sector_string = 's'+sector_string 

		for nffi,ffi in enumerate(FFI_files):
			if sector_string in ffi:
				ffi_idx = nffi 
				print('sector: ', sector_string)
				print('ffi_idx = ', ffi_idx)				

	#### call the animation function
	if return_arrays == True:
		times, flux = animate_TESS_FFI(filedir=TIC_filedir, filename=FFI_files[ffi_idx], save_animation=save_animation, return_arrays=return_arrays)
		return times, flux 

	else:
		animate_TESS_FFI(filedir=TIC_filedir, filename=FFI_files[ffi_idx], save_animation=save_animation, normalize_flux=normalize_flux, return_arrays=return_arrays)




def plot_corner(self, fitter='emcee', modelcode='batman', burnin_pct=0.1):
	print('calling _mp_visuals.py/plot_corner().')
	try:
		import corner
	except:
		print('could not import corner.')

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




def plot_bestmodel(self, fitter, modelcode, folded=False, burnin_pct=0.1, period=None):
	print('calling _mp_visuals.py/plot_bestmodel().')
	### THIS FUNCTION PLOTS YOUR BEST FIT LIGHT CURVE MODEL OVER THE DATA.

	if period == None:
		period = self.period 

	if folded == True:
		self.fold(period=period)

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






def fold(self, detrended='y', phase_offset=0.0, period=None):
	print('calling _mp_visuals.py/fold().')
	### this method will phase fold your light curve. 
	### first tau in the time series:

	if period == None:
		period = self.period 

	try:
		first_tau = self.taus[0]
	except:
		self.get_properties()
		first_tau = self.taus[0]

	ftidx = 0
	while first_tau < np.nanmin(np.hstack(self.times)):
		ftidx += 1
		first_tau = self.taus[ftidx]

	fold_times = ((((np.hstack(self.times) - first_tau - 0.5*period - phase_offset*period) % period) / period)) ### yields the remainder!
	fold_times = fold_times - 0.5

	if detrended == 'y':
		fold_fluxes = np.hstack(self.fluxes_detrend)
		fold_errors = np.hstack(self.errors_detrend)
	else:
		fold_fluxes = np.hstack(self.fluxes)
		fold_errors = np.hstack(self.errors)

	##### sort them
	fold_sort_idxs = np.argsort(fold_times)
	fold_times, fold_fluxes, folded_errors = fold_times[fold_sort_idxs], fold_fluxes[fold_sort_idxs], fold_errors[fold_sort_idxs]

	self.fold_times = fold_times
	self.fold_fluxes = fold_fluxes
	self.fold_errors = fold_errors





def examine_TPF(self, quarters=None, time_lims=None, detrend='y', mask_idxs=None):
	print('calling _mp_visuals.py/examine_TPF().')
	if type(quarters) == type(None):
		quarters = self.quarters 
	tpf_examiner(self.target, quarters=quarters, Tdur=self.duration_days, time_lims=time_lims, detrend=detrend, mask_idxs=mask_idxs, find_alias='y')





def genLS(self, show_plot = 'y', compute_fap='n', use_detrend='n', minP=None, maxP=None, LSquarters=None):
	print("calling _mp_visuals.py/genLS().")
	### this function generates a Lomb-Scargle Periodogram!
	LSperiods = []
	LSmaxpower_periods = []
	LSpowers = []
	LSfaps = []
	nquarters = len(self.quarters)

	if type(LSquarters) == type(None):
		LSquarters = self.quarters

	for qidx in np.arange(0,nquarters,1):
		this_quarter = self.quarters[qidx]
		if this_quarter not in LSquarters: ### use this in case you've only specified select quarters.
			continue

		print("processing LS for ", this_quarter)
		if nquarters != 1:
			if use_detrend == 'n':
				qtimes, qfluxes, qerrors = self.times[qidx], self.fluxes[qidx], self.errors[qidx]
			elif use_detrend == 'y':
				qtimes, qfluxes, qerrors = self.times[qidx], self.fluxes_detrend[qidx], self.errors_detrend[qidx]
		else:
			if use_detrend == 'n':
				qtimes, qfluxes, qerrors = self.times, self.fluxes, self.errors 
			elif use_use_detrend == 'y':
				qtimes, qfluxes, qerrors = self.times, self.fluxes_detrend, self.errors_detrend
		if maxP == None:
			maxperiod = 0.5 * (np.nanmax(qtimes) - np.nanmin(qtimes))
		else:
			maxperiod = maxP

		if minP == None:
			minperiod = 0.5
		else:
			minperiod = minP 
		minfreq, maxfreq = 1/maxperiod, 1/minperiod
		qls = LombScargle(qtimes, qfluxes, qerrors)
		qfreq, qpower = qls.autopower(minimum_frequency=minfreq, maximum_frequency=maxfreq)
		qperiods = 1/qfreq
		if compute_fap == 'y':
			qfap = qls.false_alarm_probability(qpower.max(), method='bootstrap')
			probabilities = [0.1, 0.05, 0.01]
			quarter_FALs = qls.false_alarm_level(probabilities)

		if show_plot == 'y':
			random_color = np.random.rand(3)
			plt.plot(qperiods[::-1], qpower[::-1], c=random_color)
			if compute_fap == 'y':
				plt.plot(qperiods[::-1], np.linspace(quarter_FALs[1], quarter_FALs[1], len(qperiods[::-1])), c=random_color)

		LSperiods.append(qperiods)
		max_power_period = qperiods[np.nanargmax(qpower)]
		LSmaxpower_periods.append(max_power_period)
		LSpowers.append(qpower)
		if compute_fap == 'y':
			LSfaps.append(qfap)

	if show_plot == 'y':
		plt.xscale('log')
		plt.xlabel('Period [days]')
		plt.ylabel('Power')
		plt.title(self.target)
		plt.show()

	LSperiods, LSpowers, LSfaps = np.array(LSperiods), np.array(LSpowers), np.array(LSfaps)

	print('LS max power periods = ', LSmaxpower_periods)
	LSperiod_median = np.nanmedian(LSmaxpower_periods)
	LSperiod_std = np.nanstd(LSmaxpower_periods)
	print('median(LS max power periods) = ', LSperiod_median)
	print('std(LS max power periods) = ', LSperiod_std)

	self.LSperiods = LSperiods
	self.LSpowers = LSpowers 
	self.LSfaps = LSfaps
	self.LSmaxperiods = LSmaxpower_periods 





def correlated_noise_detector(self):
	print('calling _mp_visuals.py/correlated_noise_detector().')
	#### this function will take the light curve and bin it up in larger and larger bins, to test whether there is correlated noise present.
	#### must mask transits (use self.mask_transit_idxs, shape=self.times.shape).

	unmasked_times, unmasked_fluxes = [], []

	for i in np.arange(0,self.times.shape[0],1):
		if len(self.mask_transit_idxs[i]) > 0:
			unmasked_times.append(np.delete(self.times[i], self.mask_transit_idxs[i]))
			unmasked_fluxes.append(np.delete(self.fluxes_detrend[i], self.mask_transit_idxs[i]))
		else:
			unmasked_times.append(self.times[i])
			unmasked_fluxes.append(self.fluxes_detrend[i])

	unmasked_times, unmasked_fluxes = np.array(unmasked_times), np.array(unmasked_fluxes)
	unmasked_times, unmasked_fluxes = np.concatenate(unmasked_times), np.concatenate(unmasked_fluxes)

	bin_sizes = np.logspace(1,4,100) ### from 1 to 10,000 per bin.
	bin_size_stds = []
	for bin_size in bin_sizes:
		print('bin_size = ', bin_size)
		#### number of bins will be length(unmasked_fluxes) / bin_size
		num_bins = int(len(unmasked_fluxes) / bin_size) + 1 
		bin_medians = []
		for bin_number in np.arange(0,num_bins-1,1):
			bin_vals = unmasked_fluxes[int(bin_number*bin_size):int((bin_number+1)*bin_size)]
			bin_medians.append(np.nanmedian(bin_vals))
		bin_std = np.nanstd(bin_medians)
		bin_size_stds.append(bin_std)

	plt.plot(bin_sizes, np.array(bin_size_stds)*1e6, c='r')
	#plt.plot(bin_sizes, (-1/2)*np.log10(bin_sizes) + bin_size_stds[0], c='k')
	plt.xscale('log')
	plt.yscale('log')
	plt.ylabel(r'$\log_{10} \, \sigma$ [ppm]')
	plt.show()

	self.bin_sizes, self.bin_size_stds = bin_sizes, bin_size_stds












