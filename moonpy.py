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


##### MAJOR UPDATE -- APRIL 12, 2021
###### CLASS METHODS ARE IMPORTED FROM _mp_visuals.py, _mp_attributes.py, and _mp_manipulation.py.
######## JUST A BIT CLEANER.



moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/moonpy.py')]

plt.rcParams["font.family"] = 'serif'

hostname = socket.gethostname()
if ('tethys' in hostname) and ('sinica' in hostname):
	#moonpydir = '/data/tethys/Documents/Software/MoonPy'
	central_data_dir = '/data/tethys/Documents/Central_Data/'
elif ('Alexs-Macbook') in hostname:
	#moonpydir = '/Users/hal9000/Documents/Software/MoonPy'
	central_data_dir = '/Users/hal9000/Documents/Central_Data'
elif 'umbriel' in hostname:
	#moonpydir = '/home/cal/ateachey/Documents/MoonPy'
	central_data_dir = '/home/cal/ateachey/Documents/Central_Data/'
else:
	#moonpydir = input('Please specify the MoonPy directory (or hard-code this into moonpy.py): ')
	#central_data_dir = input("Please specify a 'central data' directory (or hard-code this into moonpy.py): ")
	### store central_data within MoonPy directory
	if os.path.exists(moonpydir+'/Central_Data'):
		pass
	else:
		os.system('mkdir '+moonpydir+'/Central_Data')
		central_data_dir = moonpydir+'/Central_Data'
print('moonpydir = ', moonpydir)
print('Light curves will be stored in '+central_data_dir)


"""
This is the MoonPy master script! To open, you should only have to type 'import moonpy'
This package is designed to do the following:
1) download light curves from Kepler, K2 or TESS (kplr, k2plr, tess)
2) generate moon models based on user specifications
3) detrend data using CoFiAM or untrendy
4) fit a model to the data using MultiNest or emcee
5) visualize the results
"""

class MoonpyLC(object):
	from _mp_visuals import plot_lc, fold, plot_corner, plot_bestmodel, examine_TPF, genLS
	from _mp_attributes import find_transit_quarters, find_aliases, get_coords, find_planet_row, get_properties, find_taus, mystery_solver, find_neighbors, get_neighbors, find_TTVs
	from _mp_manipulation import fit, prep_for_CNN, initialize_priors, detrend, gen_batman


	### this is the light curve object. You can detrend this light curve, or fit it.
	### when you initialize it, you'll either give it the times, fluxes, and errors, OR
	### you'll provide a targetID and telescope, which will allow you to download the dataset!

	def __init__(self, targetID=None, target_type=None, lc_times=None, lc_fluxes=None, lc_errors=None, lc_flags=None, lc_quarters=None, usr_dict=None, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', sc=False, ffi='n', save_lc='y', load_lc='n', download='y', is_neighbor='n', attributes_only='n', clobber=None):
		
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
		#download_directory = central_data_dir+'Kepler_lightcurves/KIC'+str(query_format_number)
		#savepath = central_data_dir+'/'
		if self.telescope.lower() == 'kepler':
			savepath = kepler_URL_generator(find_KIC_alias(targetID))[2] 
			#savepath = central_data_dir+'/Kepler_lightcurves/'+targetID
		elif self.telescope.lower() == 'k2':
			#savepath = central_data_dir+'/K2_lightcurves/'+targetID
			savepath = k2_URL_generator(find_EPIC_alias(targetID))[2]
		elif self.telescope.lower() == 'tess':
			savepath = central_data_dir+'/TESS_lightcurves/'+targetID

		if os.path.exists(savepath):
			pass
		else:
			os.system('mkdir '+savepath)	

		### make it an attribute
		self.savepath = savepath




		if (load_lc == 'n') and (clobber == None):
			### check to see if a file already exists!
			if os.path.exists(self.savepath+'/'+str(targetID)+'_'+self.telescope+'_lightcurve.tsv'):
				print(self.savepath+'/'+str(targetID)+'_'+self.telescope+'_lightcurve.tsv exists.')
				clobber = input('Clobber? y/n: ')
				if clobber == 'n':
					load_lc = 'y'
				else:
					print('loading '+self.savepath+'/'+str(targetID)+'_'+self.telescope+'_lightcurve.tsv')

		elif (load_lc == 'n') and (clobber == 'n'):
			load_lc = 'y'
		elif (load_lc == 'n') and (clobber == 'y'):
			pass 

		if load_lc == 'y':
			save_lc = 'n' ### don't re-write what you've already got!

		if load_lc == 'y':
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
					pandafile = pandas.read_csv(self.savepath+'/'+target_name+'_'+self.telescope+'_lightcurve.tsv', delimiter='\t')
				
				except:
				
					try:
						pandafile = pandas.read_csv(self.savepath+'/'+target_name+'_lightcurve.tsv', delimiter='\t') ### older files lacked telescope information.
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
					except:
						print("could not load detrended fluxes.")
						pfluxes_detrend = pfluxes
						perrors_detrend = perrors 

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
		
					if len(lc_quarters) > 1:
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = np.array(lc_times, dtype=object), np.array(lc_fluxes, dtype=object), np.array(lc_errors, dtype=object), np.array(lc_flags, dtype=object), np.array(lc_quarters, dtype=object)
					else:
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags), np.array(lc_quarters)
					self.times, self.fluxes, self.errors, self.flags, self.quarters = lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters

					try:
						if len(lc_quarters) > 1:
							lc_fluxes_detrend, lc_errors_detrend = np.array(lc_fluxes_detrend, dtype=object), np.array(lc_errors_detrend, dtype=object)
						else:
							lc_fluxes_detrend, lc_errors_detrend = np.array(lc_fluxes_detrend), np.array(lc_errors_detrend)

						self.fluxes_detrend, self.errors_detrend = lc_fluxes_detrend, lc_errors_detrend 
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
				mast_rowidx = np.nan 
				mast_data = ascii.read('confirmed_planets.txt')
				mast_columns = mast_data.columns
				exofop_data = pandas.read_csv('exofop_toilists.pipe', delimiter='|')
				exofop_columns = exofop_data.columns



			### KEPLER HANDLING
			elif telescope.lower() == 'kepler':
				print('downloading via kplr...')
				mast_rowidx, mast_data, NEA_targetname = self.find_planet_row(row_known='n') ### cannot access ExoFOP for Kepler without a login.
				self.NEA_targetname = NEA_targetname
				try:
					#if (type(lc_times) == type(None)) and (type(lc_fluxes) == type(None)) and (type(lc_errors) == type(None)):
					if user_supplied == 'n':
						if download == 'y':
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, clobber='n', targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
						elif download == 'n':
							print('Assuming this is a neighbor... using the same times, fluxes, errors, flags and quarters!')
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = self.times, self.fluxes, self.errors, self.flags, self.quarters
					print('lc_quarters = ', lc_quarters)

				except:
					try:
						### maybe it needs the full name!
						#if (type(lc_times) == type(None)) and (type(lc_fluxes) == type(None)) and (type(lc_errors) == type(None)):
						if user_supplied == 'n':
							if download == 'y':
								lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, clobber='n', targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)	
							elif download == 'n':
								print('Assuming this is a neighbor... using the same times, fluxes, errors, flags and quarters!')
								lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = self.times, self.fluxes, self.errors, self.flags, self.quarters

						print('lc_quarters = ', lc_quarters)

					except:
						traceback.print_exc()
				print('downloaded.')





			### K2 HANDLING
			elif telescope.lower() == 'k2':
				#try:
				mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row(row_known='n') ### exofop data is available for K2 targets without a login.
				self.NEA_targetname = NEA_targetname

				print('mast_rowidx = ', mast_rowidx)
				print('exofop_rowidx = ', exofop_rowidx)
				print('NEA_targetname = ', NEA_targetname)

				#except:
				#traceback.print_exc()
				try:
					print('first try statement...')
					#if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
					lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, clobber='n', targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
				except:
					try: ### maybe it just wants the number.
						print('second try statement...')
						if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, clobber='n', targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
					except:
						traceback.print_exc()

				print('k2 lc_quarters = ', lc_quarters)





			### TESS HANDLING
			elif telescope.lower() == 'tess':
				mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row(row_known='n') ### exofop data is available for TESS targets without a login.
				self.NEA_targetname = NEA_targetname

				if ffi == 'y': ### eleanor is designed to download FFI light curves.
					lc_times, lc_fluxes, lc_errors = eleanor_target_download(targetID, lc_format=lc_format, sectors=quarters)

				elif ffi == 'n':
					if self.target.lower().startswith("toi"):
						ticnum = int(np.array(exofop_data['TIC ID'])[exofop_rowidx]) ### identify the TIC# based on the matching TOI row.
						ticname = 'TIC '+str(ticnum)
						print('calling tess_target_download().')
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(ticname)	
						print('lc_times.shape = ', lc_times.shape)

					elif self.target.lower().startswith('tic'):
						print('calling tess_target_download().')
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(targetID) ### TIC number -- LACKS the TIC prefix!
						print('lc_times.shape = ', lc_times.shape)

					else: ### for targets that don't have a TIC number or a TOI number, have to find the tic number using mast_data!
						### rowidx corresponds to the mast_data!

						### first, see if your planet name is in the comments!
						try:
							### does not provide a TIC number!!!
							tess_ra, tess_dec = np.array(mast_data['ra'])[mast_rowidx][0], np.array(mast_data['dec'])[mast_rowidx][0]
						except:
							tess_ra, tess_dec = np.array(exofop_data['RA'])[exofop_rowidx], np.array(exofop_data['Dec'])[exofop_rowidx]
						
						ticnum = np.array(exofop_data['TIC ID'][exofop_rowidx])
						ticname = 'TIC '+str(ticnum)

						print('tess_ra, tess_dec ', tess_ra, tess_dec)
						self.RA, self.Dec = tess_ra, tess_dec 
						try:
							print('calling tess_target_download().')
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(ticname)
						except:
							print('calling tess_coord_download().')
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, simbad_name = tess_coord_download(tess_ra, tess_dec)

						#### if lc_times is empty, you need to find its alias and query via tess_coord_download

			#self.rowidx = rowidx #### IMPORTANT: only define self.rowidx for your TARGET! not for neighbors!
			self.mast_rowidx = mast_rowidx

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
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
				
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
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters, target_name = kplr_coord_download(RA, Dec, coord_format=coord_format, quarters=quarters, search_radius=search_radius, lc_format=lc_format, sc=sc)
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
		#for qidx,quarters in enumerate(lc_quarters):
		for qidx in np.arange(0,nquarters,1):
			### remove NaNs
			#if nquarters != 1:
			try:
				nan_idxs = np.where(np.isfinite(np.array(lc_fluxes, dtype=np.float64)[qidx]) == False)[0]
			except:
				nan_idxs = np.where(np.isfinite(lc_fluxes[qidx]) == False)[0]
			#elif nquarters == 1:
			#	nan_idxs = np.where(np.isfinite(lc_fluxes) == False)[0]

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


			#if load_lc == 'n':
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
				#self.get_properties(locate_neighbor='y')
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
				lcfile = open(self.savepath+'/'+str(target_name)+'_'+self.telescope+'_lightcurve.tsv', mode='w')
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
	try:
		self.get_properties()
	except:
		print('could not get_properties().')
	try:
		self.find_transit_quarters()
	except:
		print('could not find_transit_quarters().')



