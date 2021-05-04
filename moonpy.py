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


"""
This is the MoonPy master script! To open, you should only have to type 'import moonpy'
This package is designed to do the following:
1) download light curves from Kepler, K2 or TESS (kplr, k2plr, tess)
2) generate moon models based on user specifications
3) detrend data using CoFiAM or untrendy
4) fit a model to the data using MultiNest or emcee
5) visualize the results

##### MAJOR UPDATE -- APRIL 12, 2021
###### CLASS METHODS ARE IMPORTED FROM _mp_visuals.py, _mp_attributes.py, and _mp_manipulation.py.
######## JUST A BIT CLEANER.
"""


first_kepler = 'y'
first_koi = 'y'
first_k2 = 'y'
first_tess = 'y'


moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/moonpy.py')]

plt.rcParams["font.family"] = 'serif'

hostname = socket.gethostname()
if ('tethys' in hostname) and ('sinica' in hostname):
	central_data_dir = '/data/tethys/Documents/Central_Data'
elif ('Alexs-Macbook') in hostname:
	central_data_dir = '/Users/hal9000/Documents/Central_Data'
elif 'umbriel' in hostname:
	central_data_dir = '/home/cal/ateachey/Documents/Central_Data'
else:
	### store central_data within MoonPy directory
	if os.path.exists(moonpydir+'/Central_Data'):
		pass
	else:
		os.system('mkdir '+moonpydir+'/Central_Data')
		central_data_dir = moonpydir+'/Central_Data'
print('moonpydir = ', moonpydir)
print('Light curves will be stored in '+central_data_dir)



class MoonpyLC(object):
	from _mp_visuals import plot_lc, fold, plot_corner, plot_bestmodel, examine_TPF, genLS, correlated_noise_detector
	from _mp_attributes import find_transit_quarters, find_aliases, get_coords, find_planet_row, get_properties, get_databases, find_taus, mystery_solver, find_neighbors, find_TTVs
	from _mp_manipulation import fit, prep_for_CNN, initialize_priors, detrend, gen_batman


	### this is the light curve object. You can detrend this light curve, or fit it.
	### when you initialize it, you'll either give it the times, fluxes, and errors, OR
	### you'll provide a targetID and telescope, which will allow you to download the dataset!

	def __init__(self, targetID=None, target_type=None, lc_times=None, lc_fluxes=None, lc_errors=None, lc_flags=None, lc_quarters=None, usr_dict=None, mask_multiple=5, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', short_cadence=False, ffi='n', save_lc='y', load_lc='n', download='y', is_neighbor='n', attributes_only='n', clobber=None):
		

		global first_kepler, first_k2, first_tess, first_koi

		self.mask_multiple = mask_multiple 

		### FOR A USER-GENERATED LIGHT CURVE, DO EVERYTHING UP TOP!
		### treat the times, fluxes and errors as a single quarter
		original_target_input = targetID ### keep this around! sometimes you want it!


		if attributes_only == 'y':
			clobber = 'n' #### if the file already exists, don't clobber it!
			##### this option allows you to download all the attributes of the object without getting the light curve!	
			###### to pull this off you will need to act as though these are user suppled light curves.
			print('attributes_only option selected... creating dummy light curves.')
			lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags, lc_quarters = np.linspace(0,1,1000), np.linspace(1,1,1000), np.linspace(0.001,0.001,1000), np.linspace(1,1,1000), np.linspace(0.001,0.001, 1000), np.linspace(0,0,1000), np.array([1])

		if (type(lc_times) != type(None)) and (type(lc_fluxes) != type(None)) and (type(lc_errors) != type(None)):
			print('using USER-SUPPLIED VALUES.')
			user_supplied = 'y'
			### if you've supplied times, fluxes, and errors 
			self.times, self.fluxes, self.errors = lc_times, lc_fluxes, lc_errors
			
			if type(lc_flags) == None:
				lc_flags = np.linspace(0,0,len(self.times))
			if type(lc_quarters) == None:
				lc_quarters = np.array([1])

			self.flags = lc_flags
			self.quarters = lc_quarters 

			if attributes_only == 'n':			
				targetID = 'USR-'+str(np.random.randint(low=0, high=1e4)+round(np.random.random(), 2))
				telescope='USER'
				self.telescope = telescope

				RA, Dec = 0.0, 0.0
				is_neighbor='y'

				try:
					self.period = usr_dict['period']
					self.tau0 = usr_dict['tau0']
					self.impact = usr_dict['impact']
					self.duration_hours = usr_dict['duration_hours']
					self.rprstar = usr_dict['rprstar']
					self.sma_AU = usr_dict['sma_AU']
					self.rp_rearth = usr_dict['rp_rearth']
				except:
					print(usr_dict)
					self.period = float(input('Enter the period: '))
					self.tau0 = float(input('Enter the tau0: '))
					self.impact = float(input('Enter the impact parameter: '))
					self.duration_hours = float(input('Enter the duration of the transit in hours: '))
					self.rprstar = float(input('Enter the Rp/Rstar: '))
					self.sma_AU = float(input('Enter the planet sma in AU: '))
					self.rp_rearth = float(input('Enter the planet radius in units of Earth radii: '))
		
		else:
			user_supplied = 'n'


		if targetID.lower().startswith('usr'):
			user_supplied = 'y'


		if type(targetID) != type(None):
			### targetID has been supplied
			if targetID.lower().startswith('usr'):
				try:
					self.period = usr_dict['period']
					self.tau0 = usr_dict['tau0']
					self.impact = usr_dict['impact']
					self.duration_hours = usr_dict['duration_hours']
					self.rprstar = usr_dict['rprstar']
					self.sma_AU = usr_dict['sma_AU']
					self.rp_rearth = usr_dict['rp_rearth']
				except:
					print(usr_dict)
					self.period = float(input('Enter the period: '))
					self.tau0 = float(input('Enter the tau0: '))
					self.impact = float(input('Enter the impact parameter: '))
					self.duration_hours = float(input('Enter the duration of the transit in hours: '))
					self.rprstar = float(input('Enter the Rp/Rstar: '))
					self.sma_AU = float(input('Enter the planet sma in AU: '))
					self.rp_rearth = float(input('Enter the planet radius in units of Earth radii: '))			

		if type(target_type) != type(None):
			self.target_type = target_type 

		if telescope == None: # if user hasn't specified, figure it out!
			if (str(targetID).lower().startswith("kic")) or (str(targetID).lower().startswith('kepler')) or (str(targetID).lower().startswith('koi')):
				telescope='kepler'
			elif str(targetID).lower().startswith('tic') or str(targetID).lower().startswith('toi') or str(targetID).lower().startswith("epic"):
				telescope='tess'
			elif (str(targetID).lower().startswith("k2")) or (str(targetID).lower().startswith('epic')):
				telescope='k2'
			else:
				telescope = input('Please specify the telescope: ')


		self.telescope = telescope 
		target_name = targetID ### maybe redundant but possibly helpful.
		self.target = target_name ### preserves the full name of the target, not just the number. 


		#### SPECIFY THE SAVEPATH FOR THE LIGHT CURVE -- NEW HANDLING (NOV 2020) SAVES THEM OUTSIDE THE MOONPY DIRECTORY, IN CENTRAL_DATA.
		if self.telescope.lower() == 'kepler':
			savepath = kepler_URL_generator(find_KIC_alias(targetID), short_cadence=short_cadence)[2] 
		elif self.telescope.lower() == 'k2':
			savepath = k2_URL_generator(find_EPIC_alias(targetID))[2]
		elif self.telescope.lower() == 'tess':
			savepath = central_data_dir+'/TESS_lightcurves/'+targetID

		savepath = nospaces(savepath)

		if os.path.exists(savepath):
			pass
		else:
			os.system('mkdir '+savepath)	

		### make it an attribute
		self.savepath = savepath


		if self.target.lower().startswith('kepler') or self.target.lower().startswith('koi') or self.target.lower().startswith('kic'):
			if self.target.lower().startswith('kepler'):
				if first_kepler == 'y':
					global kepler_NEA_data, kepler_NEA_columns, kepler_exofop_data, kepler_exofop_columns 
					kepler_NEA_data, kepler_NEA_columns, kepler_exofop_data, kepler_exofop_columns = self.get_databases()
					self.NEA_data = kepler_NEA_data
					self.NEA_columns = kepler_NEA_columns
					self.exofop_data = kepler_exofop_data
					self.exofop_columns = kepler_exofop_columns
					first_kepler = 'n'

				elif first_kepler == 'n':
					self.NEA_data = kepler_NEA_data
					self.NEA_columns = kepler_NEA_columns
					self.exofop_data = kepler_exofop_data
					self.exofop_columns = kepler_exofop_columns			

			elif (self.target.lower().startswith('koi')) or (self.target.lower().startswith('kic')):
				if first_koi == 'y':
					global koi_NEA_data, koi_NEA_columns, koi_exofop_data, koi_exofop_columns 
					koi_NEA_data, koi_NEA_columns, koi_exofop_data, koi_exofop_columns = self.get_databases()
					self.NEA_data = koi_NEA_data
					self.NEA_columns = koi_NEA_columns
					self.exofop_data = koi_exofop_data
					self.exofop_columns = koi_exofop_columns
					first_koi = 'n'

				elif first_koi == 'n':
					self.NEA_data = koi_NEA_data
					self.NEA_columns = koi_NEA_columns
					self.exofop_data = koi_exofop_data
					self.exofop_columns = koi_exofop_columns			



		elif self.target.lower().startswith('k2') or self.target.lower().startswith('epic'):
			if first_k2 == 'y':
				global k2_NEA_data, k2_NEA_columns, k2_exofop_data, k2_exofop_columns 
				k2_NEA_data, k2_NEA_columns, k2_exofop_data, k2_exofop_columns = self.get_databases()
				self.NEA_data = k2_NEA_data
				self.NEA_columns = k2_NEA_columns
				self.exofop_data = k2_exofop_data
				self.exofop_columns = k2_exofop_columns
				first_k2 = 'n'

			elif first_k2 == 'n':
				self.NEA_data = k2_NEA_data
				self.NEA_columns = k2_NEA_columns
				self.exofop_data = k2_exofop_data
				self.exofop_columns = k2_exofop_columns			
		

		
		elif self.target.lower().startswith('tess') or self.target.lower().startswith('tic') or self.target.lower().startswith('toi'):
			if first_tess == 'y':
				global tess_NEA_data, tess_NEA_columns, tess_exofop_data, tess_exofop_columns 
				tess_NEA_data, tess_NEA_columns, tess_exofop_data, tess_exofop_columns = self.get_databases()
				self.NEA_data = tess_NEA_data
				self.NEA_columns = tess_NEA_columns
				self.exofop_data = tess_exofop_data
				self.exofop_columns = tess_exofop_columns
				first_tess = 'n'

			elif first_tess == 'n':
				self.NEA_data = tess_NEA_data
				self.NEA_columns = tess_NEA_columns
				self.exofop_data = tess_exofop_data
				self.exofop_columns = tess_exofop_columns			



		### check to see if a file already exists!

		if (clobber == None) and (os.path.exists(self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv')):
			print(self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv exists.')
			clobber = input('Clobber? y/n: ')

		if (clobber == 'n') and (os.path.exists(self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv')):
			print(self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv exists.')
			load_lc = 'y'
			print('loading '+self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv')			
			
		elif (clobber == 'y') and (os.path.exists(self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv')):
			print('CLOBBERING '+self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv.')			
			load_lc = 'n'

		elif (clobber == 'n') and (os.path.exists(self.savepath+'/'+nospaces(str(targetID).lower())+'_'+self.telescope+'_lightcurve.tsv') == False):
			load_lc = 'n'

		else:
			load_lc = 'n'

		if load_lc == 'y':
			save_lc = 'n' ### don't re-write what you've already got!
			self.newlc = 'n'
		
		else:
			self.newlc = 'y'

		### strip off the prefix, find just the numbers and letters.
		if str(targetID).lower().startswith('kic'):
			targetID = targetID[3:]
		elif str(targetID).lower().startswith('tic'):
			targetID = targetID[3:]
		elif str(targetID).lower().startswith('kepler-'):
			targetID = targetID[7:]
		elif str(targetID).lower().startswith('koi-'):
			targetID = targetID[4:]
		elif str(targetID).lower().startswith('toi-'):
			targetID = targetID[4:]
		elif str(targetID).lower().startswith('k2-'):
			targetID = targetID[3:]
		elif str(targetID).lower().startswith("epic"):
			targetID = targetID[4:]
		elif str(targetID).lower().startswith('usr'):
			pass

		if str(targetID).startswith(' ') or str(targetID).startswith('-'):
			### isolate just the number!
			targetID = targetID[1:]

		self.targetID = targetID


		### intuit whether the targetID is a 'planet' (includes a letter), a KOI (includes a decimal), or a KIC (neither).
		if target_type==None: ### not specified.
			if '.' in str(self.target) and ((telescope.lower() =='kepler')):
				target_type='koi' ### since numbers are something like KOI-101.01
			elif '.' in str(self.target) and ((telescope.lower() =='tess')):
				target_type='toi' ### since TOIs are something like TOI-101.01
			elif (('b' in str(self.target)) or ('c' in str(self.target)) or ('d' in str(self.target)) or ('e' in str(self.target)) or ('f' in str(self.target)) or ('g' in str(self.target))) and ((telescope=='kepler') or (telescope=='Kepler')):
				target_type='planet' ### confirmed planets are given letters to identify them, a la Kepler-1625b.
			elif (('b' in str(self.target)) or ('c' in str(self.target)) or ('d' in str(self.target)) or ('e' in str(self.target)) or ('f' in str(self.target)) or ('g' in str(self.target))) and ((telescope=='k2') or (telescope=='K2')):		
				target_type='planet' ### k2 confirmed planets have names like K2-17b
			elif (('epic' in self.target.lower())) and ((telescope.lower() =='k2')):
				target_type = 'epic' 
			elif telescope.lower() == 'user':
				target_type = 'usr'
			else:
				### when in doubt, try these.
				if telescope.lower() == 'kepler':
					target_type='kic'
				elif telescope.lower() == 'k2':
					target_type = 'epic'
				elif telescope.lower() == 'tess':
					target_type = 'tic'

		self.target_type = target_type

		print("targetID = ", targetID)
		print('target_type = ', target_type)
		print('telescope = ', telescope)


		#### LOADING A LIGHT CURVE THAT'S ALREADY BEEN DOWNLOADED.
		if load_lc == 'y':
			if self.target.lower().startswith('usr'):
				self.telescope='user'
			elif self.target.lower().startswith('k2') or (self.target.lower().startswith('epic')):
				self.telescope = "k2"
			elif self.target.lower().startswith('kepler') or self.target.lower().startswith('kic') or self.target.lower().startswith('koi'):
				self.telescope = "kepler"
			elif (self.target.lower().startswith('tic')) or (self.target.lower().startswith('toi')):
				self.telescope = 'tess'
			else:
				try:
					print(telescope)
				except:
					telescope = input('Please specify the telescope: ')
				self.telescope = telescope 

			try:
				try:
					pandafile = pandas.read_csv(self.savepath+'/'+nospaces(target_name.lower())+'_'+self.telescope+'_lightcurve.tsv', delimiter='\t')
				
				except:
				
					try:
						pandafile = pandas.read_csv(self.savepath+'/'+nospaces(target_name.lower())+'_lightcurve.tsv', delimiter='\t') ### older files lacked telescope information.
					except:
						print("could not load the light curve from file. Will download.")
						load_lc = 'n'
				
				try:
					if self.telescope.lower() == 'user':
						ptimes = np.array(pandafile['BJD'])
					elif self.telescope.lower() == 'kepler' or self.telescope.lower() == 'k2':
						ptimes = np.array(pandafile['BKJD'])
					elif self.telescope.lower() == 'tess':
						
						try:
							ptimes = np.array(pandafile['BTJD'])	
						except:
							ptimes = np.array(pandafile['BKJD']) ### fix for earlier mislabeling of TESS time offsets.
					
					pfluxes = np.array(pandafile['fluxes'])
					perrors = np.array(pandafile['errors'])
					pquarters = np.array(pandafile['quarter'])
					pflags = np.array(pandafile['flags'])

					try:
						pfluxes_detrend = np.array(pandafile['fluxes_detrended'])
						perrors_detrend = np.array(pandafile['errors_detrended'])
						pdetrend_model = np.array(pandafile['detrend_model'])						
					
					except:
						print("could not load detrended fluxes.")
						pfluxes_detrend = np.linspace(np.nan, np.nan, len(pfluxes))
						perrors_detrend = np.linspace(np.nan, np.nan, len(pfluxes))
						pdetrend_model = np.linspace(np.nan, np.nan, len(pfluxes))

					unique_quarters = np.unique(pquarters)
					lc_times, lc_fluxes, lc_errors, lc_detrend_model, lc_fluxes_detrend, lc_errors_detrend, lc_flags, lc_quarters = [], [], [], [], [], [], [], []

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
							lc_detrend_model.append(pdetrend_model[uqidxs])	

						except:
							#pass
							lc_fluxes_detrend.append(np.nan)
							lc_errors_detrend.append(np.nan)
							lc_detrend_model.append(np.nan)
		
					if len(lc_quarters) > 1:
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = np.array(lc_times, dtype=object), np.array(lc_fluxes, dtype=object), np.array(lc_errors, dtype=object), np.array(lc_flags, dtype=object), np.array(lc_quarters, dtype=object)
					else:
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags), np.array(lc_quarters)
					
					self.times, self.fluxes, self.errors, self.flags, self.quarters = lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters

					try:
						if len(lc_quarters) > 1:
							lc_fluxes_detrend, lc_errors_detrend = np.array(lc_fluxes_detrend, dtype=object), np.array(lc_errors_detrend, dtype=object)
							lc_detrend_model = np.array(lc_detrend_model, dtype=object)
						else:
							lc_fluxes_detrend, lc_errors_detrend = np.array(lc_fluxes_detrend), np.array(lc_errors_detrend)
							lc_detrend_model = np.array(lc_detrend_model)

						self.fluxes_detrend, self.errors_detrend = lc_fluxes_detrend, lc_errors_detrend 
						self.detrend_model = lc_detrend_model 
					
					except:
						pass

				except:
					print("could not load the light curve from file. Will download.")
					load_lc = 'n'

			except:
				traceback.print_exc()
				print("could not load the light curve from file. Will download.")
				load_lc = 'n'



		### HANDLING FOR DOWNLOADING A FRESH LIGHT CURVE.
		if (load_lc=='n') and (type(targetID) != type(None)) and (type(telescope) != type(None)) and (attributes_only == 'n'):
			### implies you've selected a target you want to download.

			### USER-INPUT HANDLING
			if telescope.lower() == 'user':
				lc_times, lc_fluxes, lc_errors, lc_flags, lcquarters = self.times, self.fluxes, self.errors, self.flags, self.quarters
				NEA_rowidx = np.nan 
				NEA_data = ascii.read('confirmed_planets.txt')
				NEA_columns = NEA_data.columns
				exofop_data = pandas.read_csv('exofop_toilists.pipe', delimiter='|')
				exofop_columns = exofop_data.columns



			### KEPLER HANDLING
			elif telescope.lower() == 'kepler':
				print('downloading via kplr...')
				NEA_rowidx, NEA_targetname = self.find_planet_row(row_known='n') ### cannot access ExoFOP for Kepler without a login.
				self.NEA_targetname = NEA_targetname
				try:
					if user_supplied == 'n':
						if download == 'y':
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, clobber=clobber, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, short_cadence=short_cadence)
						elif download == 'n':
							print('Assuming this is a neighbor... using the same times, fluxes, errors, flags and quarters!')
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = self.times, self.fluxes, self.errors, self.flags, self.quarters
					print('lc_quarters = ', lc_quarters)


				except:
					try:
						### maybe it needs the full name!
						if user_supplied == 'n':
							if download == 'y':
								lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, clobber=clobber, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, short_cadence=short_cadence)	
							elif download == 'n':
								print('Assuming this is a neighbor... using the same times, fluxes, errors, flags and quarters!')
								lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = self.times, self.fluxes, self.errors, self.flags, self.quarters

						print('lc_quarters = ', lc_quarters)

					except:
						traceback.print_exc()
				print('downloaded.')

				if np.isfinite(NEA_rowidx) == False:
					try:
						self.find_aliases()
						print("LOOKING FOR ALIASES IN NEA data.")
						for alias in self.aliases:
							alias_NEA_rowidx = self.find_planet_row(alias=alias, row_known='n')[0] ### exofop data is available for TESS targets without a login.			
							try:
								NEA_rowidx = int(alias_NEA_rowidx)
								break
							except:
								continue
					except:
						pass

				print('[mp] NEA_rowidx = ', NEA_rowidx)
				self.NEA_rowidx = NEA_rowidx 






			### K2 HANDLING
			elif telescope.lower() == 'k2':
				#try:
				NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='n') ### exofop data is available for K2 targets without a login.
				self.NEA_targetname = NEA_targetname

				print('NEA_rowidx = ', NEA_rowidx)
				print('exofop_rowidx = ', exofop_rowidx)
				print('NEA_targetname = ', NEA_targetname)

				try:
					print('first try statement...')
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, clobber=clobber, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, short_cadence=short_cadence)
				except:
					try: ### maybe it just wants the number.
						print('second try statement...')
						if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, clobber=clobber, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, short_cadence=short_cadence)
					except:
						traceback.print_exc()

				print('k2 lc_quarters = ', lc_quarters)

				if np.isfinite(NEA_rowidx) == False:
					try:
						self.find_aliases()
						print("LOOKING FOR ALIASES IN NEA data.")
						for alias in self.aliases:
							alias_NEA_rowidx = self.find_planet_row(alias=alias, row_known='n')[0] ### exofop data is available for TESS targets without a login.			
							try:
								NEA_rowidx = int(alias_NEA_rowidx)
								break
							except:
								continue
					except:
						pass


				if np.isfinite(exofop_rowidx) == False:
					try:
						self.find_aliases()
						print('LOOKING FOR ALIASES in ExoFOP data.')
						for alias in self.aliases:
							alias_exofop_rowidx = self.find_planet_row(alias=alias, row_known='n')[1] ### exofop data is available for TESS targets without a login.			
							try:
								exofop_rowidx = int(alias_exofop_rowidx)
								break
							except:
								continue
					except:
						pass


				print('[mp] NEA_rowidx = ', NEA_rowidx)
				print('[mp] exofop_rowidx = ', exofop_rowidx)
				self.NEA_rowidx = NEA_rowidx
				self.exofop_rowidx = exofop_rowidx






			### TESS HANDLING
			elif telescope.lower() == 'tess':
				try:
					NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='n') ### exofop data is available for TESS targets without a login.			
				except:
					print('UNABLE TO FIND NEA_rowidx and/or exofop_rowidx. May not exist.')
					NEA_rowidx = np.nan 
					exofop_rowidx = np.nan 
					NEA_targetname = targetID

				try:
					NEA_rowidx = int(NEA_rowidx)
				except:
					pass

				try:
					exofop_rowidx = int(exofop_rowidx)
				except:
					pass

				if np.isfinite(NEA_rowidx) == False:
					try:
						self.find_aliases()
						print("LOOKING FOR ALIASES IN NEA data.")
						for alias in self.aliases:
							alias_NEA_rowidx = self.find_planet_row(alias=alias, row_known='n')[0] ### exofop data is available for TESS targets without a login.			
							try:
								NEA_rowidx = int(alias_NEA_rowidx)
								break
							except:
								continue
					except:
						pass

				if np.isfinite(exofop_rowidx) == False:
					try:
						self.find_aliases()
						print('LOOKING FOR ALIASES in ExoFOP data.')
						for alias in self.aliases:
							alias_exofop_rowidx = self.find_planet_row(alias=alias, row_known='n')[1] ### exofop data is available for TESS targets without a login.			
							try:
								exofop_rowidx = int(alias_exofop_rowidx)
								break
							except:
								continue
					except:
						pass


				print('[mp] NEA_rowidx = ', NEA_rowidx)
				print('[mp] exofop_rowidx = ', exofop_rowidx)
				self.NEA_rowidx = NEA_rowidx
				self.exofop_rowidx = exofop_rowidx
				self.NEA_targetname = NEA_targetname



				if ffi == 'n':

					if self.target.lower().startswith("toi"):
						try:
							ticnum = int(np.array(exofop_data['TIC ID'])[self.exofop_rowidx]) ### identify the TIC# based on the matching TOI row.
							ticname = 'TIC '+str(ticnum)
						
						except:
							try:
								ticname = find_TIC_alias(self.target)
								ticnum = ticname[ticname.find('TIC'):]
								if ticnum.startswith(' ') or ticnum.startswith('-'):
									ticnum = ticnum[1:]
							except:
								print('Could not access TIC number in NEA catalog.')

						try:
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(ticname)	
							print('lc_times.shape = ', lc_times.shape)
						except:
							print('Unable to download lc with tess_target_download(). Likely missing TIC.')



					elif self.target.lower().startswith('tic'):
						print('calling tess_target_download().')
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(targetID) ### TIC number -- LACKS the TIC prefix!
						print('lc_times.shape = ', lc_times.shape)



					else: 
						### for targets that don't have a TIC number or a TOI number, have to find the tic number using NEA_data!
						### rowidx corresponds to the NEA_data!
						### first, see if your planet name is in the comments!

						try:
							### does not provide a TIC number!!!
							tess_ra, tess_dec = np.array(NEA_data['ra'])[NEA_rowidx][0], np.array(NEA_data['dec'])[NEA_rowidx][0]
							ticname = find_TIC_alias(self.target)
							ticnum = ticname[ticname.find('TIC'):]
							if ticnum.startswith(' ') or ticnum.startswith('-'):
								ticnum = ticnum[1:]							
						except:
							try:
								tess_ra, tess_dec = np.array(NEA_data['ra'])[NEA_rowidx], np.array(NEA_data['dec'])[NEA_rowidx]
								ticname = find_TIC_alias(self.target)
								ticnum = ticname[ticname.find('TIC'):]
								if ticnum.startswith(' ') or ticnum.startswith('-'):
									ticnum = ticnum[1:]								
							except:
								try:
									tess_ra, tess_dec = np.array(exofop_data['RA'])[exofop_rowidx], np.array(exofop_data['Dec'])[exofop_rowidx]
									ticnum = np.array(exofop_data['TIC ID'][exofop_rowidx])
									ticname = 'TIC '+str(ticnum)							

								except:								
									ticname = find_TIC_alias(self.target)
									ticnum = ticname[ticname.find('TIC'):]
									if ticnum.startswith(' ') or ticnum.startswith('-'):
										ticnum = ticnum[1:]

						try:
							print('tess_ra, tess_dec ', tess_ra, tess_dec)
							self.RA, self.Dec = tess_ra, tess_dec 
						except:
							pass

						try:
							print("ticname: ", ticname)
						except:
							pass

						try:
							print('calling tess_target_download().')
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(ticname)
						except:
							try:
								print('calling tess_coord_download().')
								lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, simbad_name = tess_coord_download(tess_ra, tess_dec)
							except:
								ffi = 'y' #### try to find it in the FULL FRAME IMAGES!


				if ffi == 'y': ### eleanor is designed to download FFI light curves.
					#lc_times, lc_fluxes, lc_errors = eleanor_target_download(targetID, lc_format=lc_format, sectors=quarters)
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = TESS_QLP_load(targetID, sectors=quarters, clobber=clobber)


				self.times, self.fluxes, self.errors, self.flags, self.quarters = lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters		


			if (self.telescope.lower() != 'kepler') and (self.telescope.lower() != 'user'):
				self.exofop_rowidx = exofop_rowidx
			self.telescope = telescope
			### SIMBAD_QUERY
			self.get_coords()

			try:
				self.find_aliases()
			except:
				print('ALIASES COULD NOT BE FOUND.')


		### HANDLING FOR USER-SUPPLIED COORDINATES.
		elif (load_lc == 'n') and (type(RA) != type(None)) and (type(Dec) != type(None)) and (attributes_only == 'n'): 
			try:	
				if (telescope.lower() == 'kepler') or (telescope.lower() =='k2'):
					if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, short_cadence=short_cadence)
				
				elif telescope == 'tess':
					print('calling tess_coord_download().')
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = tess_coord_download(RA, Dec)

				elif telescope.lower() == 'user':
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = self.times, self.fluxes, self.errors, self.flags, self.quarters
					target_name = targetID 
			
			except:
				telescope = input('Specify the telescope: ')
				self.telescope = telescope 
				if (telescope.lower() == 'kepler') or (telescope.lower() =='k2'):
					if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, short_cadence=short_cadence)
				elif telescope == 'tess':
					try:
						print('calling tess_coord_download().')
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_coord_download(RA, Dec)
					except:
						lc_times, lc_fluxes, lc_errors = eleanor_coord_download(RA, Dec)
				elif telescope.lower() == 'user':
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = self.times, self.fluxes, self.errors, self.flags, self.quarters
					target_name = targetID 
					
			self.target = target_name


		### YOU HAVEN'T DOWNLOADED A LIGHT CURVE OR SUPPLIED ONE, SO WHAT ARE YOU DOING?
		elif load_lc == 'y':
			pass

		else:
			if attributes_only == 'n':
				raise Exception("You have supplied inconsistent inputs. Must be 1) lc_times, \
					lc_fluxes, and lc_errors, 2) targetID and telescope, or 3) RA, Dec and telescope.")


		print ('load_lc = ', load_lc)



		### BELOW THIS LINE IS HANDLING FOR EVERY DIFFERENT WAY YOU MIGHT HAVE LOADED THE LIGHT CURVE ABOVE.
		if load_lc == 'y':
			pass ### you've already turned them into arrays.
		else:
			lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.array(lc_times, dtype=object), np.array(lc_fluxes, dtype=object), np.array(lc_errors, dtype=object), np.array(lc_fluxes, dtype=object), np.array(lc_errors, dtype=object), np.array(lc_flags, dtype=object)
		
		nquarters = len(lc_quarters)
		for qidx in np.arange(0,nquarters,1):
			### remove NaNs
			try:
				nan_idxs = np.where(np.isfinite(np.array(lc_fluxes, dtype=np.float64)[qidx]) == False)[0]
			except:
				nan_idxs = np.where(np.isfinite(lc_fluxes[qidx]) == False)[0]

			if len(nan_idxs) > 0:
				if len(lc_quarters) > 1:
					lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], nan_idxs), np.delete(lc_fluxes[qidx], nan_idxs), np.delete(lc_errors[qidx], nan_idxs), np.delete(lc_fluxes_detrend[qidx], nan_idxs), np.delete(lc_errors_detrend[qidx], nan_idxs), np.delete(lc_flags[qidx], nan_idxs)
					assert np.all(np.isfinite(lc_fluxes[qidx]))
					assert np.all(np.isfinite(lc_errors[qidx]))

				elif len(lc_quarters) == 1:
					lc_times = lc_times[qidx]
					lc_fluxes = lc_fluxes[qidx]
					lc_errors = lc_errors[qidx]
					lc_fluxes_detrend = lc_fluxes_detrend[qidx]
					lc_errors_detrend = lc_errors_detrend[qidx]
					lc_flags = lc_flags[qidx]
					lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.delete(lc_times, nan_idxs), np.delete(lc_fluxes, nan_idxs), np.delete(lc_errors, nan_idxs), np.delete(lc_fluxes_detrend, nan_idxs), np.delete(lc_errors_detrend, nan_idxs), np.delete(lc_flags, nan_idxs)
					assert np.all(np.isfinite(lc_fluxes))
					assert np.all(np.isfinite(lc_errors))

			if remove_flagged == 'y':
				if len(lc_quarters) > 1:
					flag_idxs = np.where(lc_flags[qidx] != 0)[0]
				elif len(lc_quarters) == 1:
					flag_idxs = np.where(lc_flags != 0)[0]

				if len(flag_idxs) > 0:
					if len(lc_quarters) > 1:
						lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = np.delete(lc_times[qidx], flag_idxs), np.delete(lc_fluxes[qidx], flag_idxs), np.delete(lc_errors[qidx], flag_idxs), np.delete(lc_fluxes_detrend[qidx], flag_idxs), np.delete(lc_errors_detrend[qidx], flag_idxs), np.delete(lc_flags[qidx], flag_idxs)
					elif len(lc_quarters) == 1:
						lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.delete(lc_times, flag_idxs), np.delete(lc_fluxes, flag_idxs), np.delete(lc_errors, flag_idxs), np.delete(lc_fluxes_detrend, flag_idxs), np.delete(lc_errors_detrend, flag_idxs), np.delete(lc_flags, flag_idxs)					


			if len(lc_quarters) > 1:
				timesort = np.argsort(lc_times[qidx])

				if lc_fluxes_detrend.shape == 0:
					lc_fluxes_detrend = lc_fluxes
				if lc_errors_detrend.shape == 0:
					lc_errors_detrend = lc_errors 
				lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = lc_times[qidx][timesort], lc_fluxes[qidx][timesort], lc_errors[qidx][timesort], lc_fluxes_detrend[qidx][timesort], lc_errors_detrend[qidx][timesort], lc_flags[qidx][timesort]

			elif len(lc_quarters) == 1:
				try:
					timesort = np.argsort(np.array(lc_times))
					lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.array(lc_times)[timesort], np.array(lc_fluxes)[timesort], np.array(lc_errors)[timesort], np.array(lc_fluxes_detrend)[timesort], np.array(lc_errors_detrend)[timesort], np.array(lc_flags)[timesort]
				except:
					try: ### in some cases the data is loaded as nested lists.
						timesort = timesort[0]
						lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = lc_times[0], lc_fluxes[0], lc_errors[0], lc_fluxes_detrend[0], lc_errors_detrend[0], lc_flags[0]
						lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = lc_times[timesort], lc_fluxes[timesort], lc_errors[timesort], lc_fluxes_detrend[timesort], lc_errors_detrend[timesort], lc_flags[timesort]
					except:
						traceback.print_exc()
						raise Exception('Something wrong with the way the data has been loaded.')
						


		self.times = lc_times
		self.fluxes = lc_fluxes
		self.errors = lc_errors
		self.flags = lc_flags
		self.quarters = lc_quarters

		### calculate which quarters have transits!
		try:
			self.find_transit_quarters()
			if is_neighbor == 'y':
				self.get_properties(locate_neighbor='n')
				self.find_transit_quarters(locate_neighbor='n') ### get_properties() is called first thing!
				traceback.print_exc()
			
			elif is_neighbor == 'n':
				self.find_transit_quarters(locate_neighbor='y') ### ditto above.
				traceback.print_exc()
				self.get_neighbors(save_to_file='n') ### do it once when you initialize, so you don't have to do it again!
		except:
			traceback.print_exc()
			print(" ")
			print(" ")
			print('ONE OR MORE METHODS FAILED. Identifier is likely not in SIMBAD, and aliases are unknown.')
			print(" ")
			print("Try calling self.mystery_solver(tau0, period, duration_hours) to manually input critical variables.")
			print("Altenatively, try the input again with a different target ID.")

		try:
			if save_lc == 'y':
				### write to a file!
				lcfile = open(self.savepath+'/'+nospaces(str(target_name).lower())+'_'+self.telescope+'_lightcurve.tsv', mode='w')
				if self.telescope.lower() == 'user':
					lcfile.write('BJD\tfluxes\terrors\tflags\tquarter\n')
				elif self.telescope.lower() == 'kepler' or self.telescope.lower() == 'k2':
					lcfile.write('BKJD\tfluxes\terrors\tflags\tquarter\n')
				elif self.telescope.lower() == 'tess':
					lcfile.write('BTJD\tfluxes\terrors\tflags\tquarter\n')				
				if len(self.quarters) > 1:
					for qidx in np.arange(0,len(self.quarters),1):
						qtq = lc_quarters[qidx]
						qtimes, qfluxes, qerrors, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_flags[qidx]
						for qt, qf, qe, qfl in zip(qtimes, qfluxes, qerrors, qflags):
							lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qfl)+'\t'+str(qtq)+'\n')
				elif len(self.quarters) == 1:
					if self.times.shape[0] == 1:
						lc_times, lc_fluxes, lc_errors = lc_times[0], lc_fluxes[0], lc_errors[0]
					else:
						pass

					for qt, qf, qe, qfl in zip(lc_times, lc_fluxes, lc_errors, lc_flags):
						lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qfl)+'\t'+str(self.quarters[0])+'\n')
				lcfile.close()
		except:
			traceback.print_exc()
			raise Exception('an exception was raised while trying to save the light curve.')



	### try this at the end
	#try:
	#	self.get_properties()
	#except:
	#	print('could not get_properties().')
	#try:
	#	self.find_transit_quarters()
	#except:
	#	print('could not find_transit_quarters().')



	def get_neighbors(self, save_to_file='y', mask_multiple=None):
		print('calling moonpy.py/get_neighbors().')

		neighbor_dict = {}
		for neighbor in self.neighbors:

			###
			try:
				if (self.telescope.lower() == 'kepler'):
				
					if neighbor.lower().startswith('kepler'):
						neighbor_key = 'K'+str(neighbor[7:])
					else:
						### for now this is just a KOI!
						neighbor_key = neighbor
						while neighbor_key.lower().startswith('k') or neighbor_key.startswith('0'):
							neighbor_key = neighbor_key[1:]
						neighbor_key = 'K'+str(neighbor_key)
				
				elif (self.telescope.lower() == 'k2'):
					neighbor_key = 'K2_'+str(neighbor_key[3:])

				elif (self.telescope.lower() == 'tess'):
					neighbor_key = neighbor 

				### now generate the LC object -- note that this will download duplicate light curves!
				### for now, just delete them after downloading... should come up with a better solution.
				print('calling MoonpyLC for neighbor: ', neighbor)

				neighbor_dict[neighbor_key] = MoonpyLC(targetID=neighbor, is_neighbor='y', clobber='n')

			except:
				traceback.print_exc()
				print("COULD NOT DOWNLOAD INFORMATION FOR "+str(neighbor))

			self.neighbor_dict = neighbor_dict 


		final_neighbor_IDs = []		
		neighbor_transit_idxs = []
		neighbor_transit_IDs = []
		
		for neighbor in neighbor_dict.keys():			
			neighbor_taus = neighbor_dict[neighbor].taus 
			neighbor_dur = neighbor_dict[neighbor].duration_days 

			print('appending '+str(len(neighbor_taus))+' neighbor transit times.')
			print('neighbor duration = '+str(neighbor_dur))

			for nt in neighbor_taus:
				ntidxs = np.where((np.hstack(self.times) >= (nt - (self.mask_multiple/2)*neighbor_dur)) & (np.hstack(self.times) <= (nt + (self.mask_multiple/2)*neighbor_dur)))[0]
				for ntidx in ntidxs:
					neighbor_transit_idxs.append(ntidx)
					neighbor_transit_IDs.append(neighbor)

		### include the target here!
		try:
			target_taus = self.taus
		except:
			print('self.taus not available. Setting target_taus = np.array([np.nan])')
			target_taus = np.array([np.nan])
		try:
			target_dur = self.duration_days
		except:
			print('self.duration_days not available. Setting target_dur = np.nan')
			target_dur = np.nan 
		for tt in target_taus:
			ttidxs = np.where((np.hstack(self.times) >= (tt - (self.mask_multiple/2)*target_dur)) & (np.hstack(self.times) <= (tt + (self.mask_multiple/2)*target_dur)))[0]
			for ttidx in ttidxs:
				neighbor_transit_idxs.append(ttidx)
				neighbor_transit_IDs.append(self.target)

		neighbor_transit_idxs = np.array(neighbor_transit_idxs)
		neighbor_transit_IDs = np.array(neighbor_transit_IDs)

		print('number of neighbor transit idxs = ', len(neighbor_transit_idxs))
		for ntidx in np.unique(neighbor_transit_idxs):
			### find all transiting planets that have this index
			all_neighbors = neighbor_transit_IDs[np.where(neighbor_transit_idxs == ntidx)[0]]
			final_neighbor_IDs.append(list(all_neighbors))


		neighbor_transit_idxs = np.unique(neighbor_transit_idxs)
		print('number of unique neighbor_transit_idxs = ', len(neighbor_transit_idxs))

		neighbor_transit_times = []
		neighbor_transit_list = []
		neighbor_transit_ID = []

		if save_to_file == 'y':
			try:
				lcfile = open(self.savepath+'/'+nospaces(self.target.lower())+"_"+self.telescope+"_lightcurve.tsv", mode='r')
			except:
				lcfile = open(self.savepath+'/'+nospaces(self.target.lower())+'_lightcurve.tsv', mode='r') ### for older files.
			lcfile_new = open(self.savepath+"/"+nospaces(self.target.lower())+"_"+self.telescope+"_lc_temp.tsv", mode='w')

			for nline, line in enumerate(lcfile):
				if nline == 0:
					newline = line[:-1]+'\tin_transit\ttransiter\n'

				else:
					if (nline-1) in neighbor_transit_idxs:
						ntidx = np.where(neighbor_transit_idxs == (nline-1))[0][0]
						every_neighbor = final_neighbor_IDs[ntidx]
						newline = line[:-1]+'\ty\t'+str(every_neighbor)+'\n'
						neighbor_transit_times.append(np.hstack(self.times)[nline-1])
						neighbor_transit_list.append('y')
						neighbor_transit_ID.append(str(every_neighbor))
					else:
						newline = line[:-1]+'\tn\n'
						neighbor_transit_list.append('n')
						neighbor_transit_ID.append('')

				lcfile_new.write(newline)
			lcfile.close()
			lcfile_new.close()
			### rename the file.
			os.system('mv '+self.savepath+'/'+nospaces(self.target.lower())+'_'+self.telescope+'_lc_temp.tsv '+self.savepath+'/'+nospaces(self.target.lower())+'_'+self.telescope+'_lightcurve.tsv')

		self.neighbor_transit_times = np.array(neighbor_transit_times)
		self.neighbor_transit_list = neighbor_transit_list
		self.neighbor_transit_IDs = neighbor_transit_ID 

