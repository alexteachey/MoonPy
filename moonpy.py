from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
#import astropy
import warnings
from astropy.io import ascii
import time
#import datetime
import pandas
import traceback
from astroquery.simbad import Simbad 
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
from astropy.stats import LombScargle

#### BELOW ARE MOONPY PACKAGES
from mp_tools import *
from mp_lcfind import *
from mp_detrend import untrendy_detrend, cofiam_detrend, george_detrend, medfilt_detrend
from mp_batman import run_batman
from mp_fit import mp_multinest, mp_emcee
from cofiam import max_order
#from pyluna import run_LUNA, prepare_files
from pyluna import prepare_files
from mp_tpf_examiner import *


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

	def __init__(self, targetID=None, target_type=None, lc_times=None, lc_fluxes=None, lc_errors=None, lc_flags=None, lc_quarters=None, usr_dict=None, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', sc=False, ffi='n', save_lc='y', load_lc='n', is_neighbor='n', clobber=None):
		### FOR A USER-GENERATED LIGHT CURVE, DO EVERYTHING UP TOP!
		### treat the times, fluxes and errors as a single quarter
		if (type(lc_times) != type(None)) and (type(lc_fluxes) != type(None)) and (type(lc_errors) != type(None)):
			print('using USER-SUPPLIED VALUES.')
			### if you've supplied times, fluxes, and errors 
			self.times, self.fluxes, self.errors = lc_times, lc_fluxes, lc_errors
			
			if type(lc_flags) == None:
				lc_flags = np.linspace(0,0,len(self.times))
			if type(lc_quarters) == None:
				lc_quarters = np.array([1])

			self.flags = lc_flags
			self.quarters = lc_quarters 

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


		if (load_lc == 'n') and (clobber == None):
			### check to see if a file already exists!
			if os.path.exists(savepath+'/'+str(targetID)+'_'+self.telescope+'_lightcurve.tsv'):
				clobber = input('light curve already exists. Clobber? y/n: ')
				if clobber == 'n':
					load_lc = 'y'

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
					pandafile = pandas.read_csv(savepath+'/'+target_name+'_'+self.telescope+'_lightcurve.tsv', delimiter='\t')
				except:
					try:
						pandafile = pandas.read_csv(savepath+'/'+target_name+'_lightcurve.tsv', delimiter='\t') ### older files lacked telescope information.
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

			except:
				traceback.print_exc()
				print("could not load the light curve from file. Will download.")
				load_lc = 'n'



		### HANDLING FOR DOWNLOADING A FRESH LIGHT CURVE.
		if (load_lc=='n') and (type(targetID) != type(None)) and (type(telescope) != type(None)):
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
					if (type(lc_times) == type(None)) and (type(lc_fluxes) == type(None)) and (type(lc_errors) == type(None)):
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
				except:
					try:
						### maybe it needs the full name!
						if (type(lc_times) == type(None)) and (type(lc_fluxes) == type(None)) and (type(lc_errors) == type(None)):
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)	
					except:
						traceback.print_exc()
				print('downloaded.')

			### K2 HANDLING
			elif telescope.lower() == 'k2':
				#try:
				mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row(row_known='n') ### exofop data is available for K2 targets without a login.
				self.NEA_targetname = NEA_targetname
				#except:
				#traceback.print_exc()
				try:
					if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(self.target, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
				except:
					try: ### maybe it just wants the number.
						if (type(lc_times) == None) and (type(lc_fluxes) == None) and (type(lc_errors) == None):
							lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, targtype=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)
					except:
						traceback.print_exc()

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
		elif (load_lc == 'n') and (type(RA) != type(None)) and (type(Dec) != type(None)): 
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
			raise Exception("You have supplied inconsistent inputs. Must be 1) lc_times, \
				lc_fluxes, and lc_errors, 2) targetID and telescope, or 3) RA, Dec and telescope.")


		print ('load_lc = ', load_lc)

		### BELOW THIS LINE IS HANDLING FOR EVERY DIFFERENT WAY YOU MIGHT HAVE LOADED THE LIGHT CURVE ABOVE.
		if load_lc == 'y':
			pass ### you've already turned them into arrays.
		else:
			lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = np.array(lc_times), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_fluxes), np.array(lc_errors), np.array(lc_flags)
		
		nquarters = len(lc_quarters)
		#for qidx,quarters in enumerate(lc_quarters):
		for qidx in np.arange(0,nquarters,1):
			### remove NaNs
			#if nquarters != 1:
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


			### sort the times here!
			"""
			print('len(lc_quarters) = ', len(lc_quarters))
			print('lc_times = ', lc_times)
			print ('lc_times.shape = ', lc_times.shape)
			print('lc_fluxes.shape = ', lc_fluxes.shape)
			print('lc_errors.shape = ', lc_errors.shape)
			print('lc_fluxes_detrend.shape = ', lc_fluxes_detrend.shape)
			print('lc_errors_detrend.shape = ', lc_errors_detrend.shape)
			print('lc_flags.shape = ', lc_flags.shape)
			print(" ")
			print(" ")
			"""

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
				lcfile = open(savepath+'/'+str(target_name)+'_'+self.telescope+'_lightcurve.tsv', mode='w')
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


	### DETRENDING!

	def detrend(self, dmeth='cofiam', save_lc='y', mask_transits='y', mask_neighbors='y', skip_ntqs='n', kernel=None, max_degree=30):
		exceptions_raised = 'n'

		if mask_neighbors == 'y':
			mask_transits = 'y' 
		### optional values for method are "cofiam", "untrendy", "medfilt"
		### EACH QUARTER SHOULD BE DETRENDED INDIVIDUALLY!

		#try:
		#	self.duration_hours ### tests whether you've called get_properties yet.
		#except:
		#	self.get_properties(locate_neighbor='n')

		master_detrend, master_error_detrend, master_flags_detrend = [], [], []

		nquarters = len(self.quarters)
		for qidx in np.arange(0,nquarters,1):
			skip_quarter = 'n'
			print('quarter = ', self.quarters[qidx])
			if nquarters != 1:
				dtimes, dfluxes, derrors, dflags = self.times[qidx], self.fluxes[qidx], self.errors[qidx], self.flags[qidx]
			elif nquarters == 1:
				dtimes, dfluxes, derrors, dflags = self.times, self.fluxes, self.errors, self.flags
			print('dtimes.shape = ', dtimes.shape)
			print('dtimes.shape[0] = ', dtimes.shape[0])
			if dtimes.shape[0] == 0:
				exceptions_raised = 'y'
				fluxes_detrend, errors_detrend = dfluxes, derrors
				flags_detrend = np.linspace(2097152,2097152,len(fluxes_detrend))
				master_detrend.append(np.array(fluxes_detrend))
				master_error_detrend.append(np.array(errors_detrend))
				master_flags_detrend.append(np.array(flags_detrend))
				continue

			#if int(dtimes.shape[0]) < 10:
			#	print('not enough times to detrend!')
			#	continue ### there's nothing to detrend here!


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
				try:
					quarter_transit_taus = self.taus[((self.taus > np.nanmin(dtimes)) & (self.taus < np.nanmax(dtimes)))]
				except:
					self.find_transit_quarters()
					traceback.print_exc()
					quarter_transit_taus = self.taus[((self.taus > np.nanmin(dtimes)) & (self.taus < np.nanmax(dtimes)))]		
								
				for qtt in quarter_transit_taus:
					in_transit_idxs = np.where((dtimes >= qtt - self.duration_days) & (dtimes <= qtt + self.duration_days))[0]
					mask_transit_idxs.append(in_transit_idxs)


				### add neighbor transit times to mask_transit_idxs.
				try:
					print(self.all_transit_times)
				except:
					try:
						self.mystery_solver(self.tau0, self.period, self.duration_hours)
					except:
						pass

				if mask_neighbors == 'y':
					if len(self.neighbors) > 0:
						neighbor_transit_idxs = []
						for ntt in self.all_transit_times:
							if ntt >= np.nanmin(dtimes) and ntt <= np.nanmax(dtimes):
								neighbor_transit_idxs.append(np.where(ntt == dtimes)[0])
						mask_transit_idxs.append(neighbor_transit_idxs)

				try:
					print("transit midtimes this quarter: ", quarter_transit_taus)
					print('min, max quarter times: ', np.nanmin(dtimes), np.nanmax(dtimes))
					if len(quarter_transit_taus) > 1:
						mask_transit_idxs = np.concatenate((mask_transit_idxs))
					else:
						mask_transit_idxs = np.array(mask_transit_idxs)

					mask_transit_idxs = np.unique(mask_transit_idxs)

					print('mask_transit_idxs = ', mask_transit_idxs)	
					if len(quarter_transit_taus) > 0:
						print("transit in this quarter.")
				except:
					traceback.print_exc()
					mask_transit_idxs = np.array([])
					print('no transits in this quarter.')
					if skip_ntqs == 'y':
						### skip this quarter! there are no transits present.
						fluxes_detrend, errors_detrend, flags_detrend = dfluxes, derrors, dflags
						skip_quarter = 'y'


			elif mask_transits == 'n':
				mask_transit_idxs = None 

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
						fluxes_detrend, errors_detrend = cofiam_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs, max_degree=max_degree)
						#flags_detrend = np.linspace(0,0,len(fluxes_detrend))
						flags_detrend = dflags

					elif dmeth == 'untrendy':
						print("UNTRENDY ISN'T REALLY SUPPORTED RIGHT NOW, SORRY!")
						fluxes_detrend, errors_detrend = untrendy_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs)
						#flags_detrend = np.linspace(0,0,len(fluxes_detrend))
						flags_detrend = dflags

					elif dmeth == 'george':
						fluxes_detrend, errors_detrend = george_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs)
						#flags_detrend = np.linspace(0,0,len(fluxes_detrend))
						flags_detrend = dflags

					elif dmeth == 'medfilt':
						print("MEDIAN FILTERING HAS YET TO BE TESTED EXTENSIVELY. BEWARE!")
						fluxes_detrend, errors_detrend = medfilt_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, size=kernel, mask_idxs=mask_transit_idxs)
						#flags_detrend = np.linspace(0,0,len(fluxes_detrend))
						flags_detrend = dflags

				except:
					traceback.print_exc()
					print('Detrending failed for this quarter. All points have likely been screened.')
					fluxes_detrend, errors_detrend = dfluxes, derrors 
					flags_detrend = np.linspace(2097152,2097152,len(fluxes_detrend))
					exceptions_raised = 'y'


			### update self -- just this quarter!
			if (skip_ntqs == 'n') and (exceptions_raised == 'n'):
				assert np.all(dfluxes != fluxes_detrend)
				assert np.all(derrors != errors_detrend)


			if len(self.quarters) > 1:
				master_detrend.append(np.array(fluxes_detrend))
				master_error_detrend.append(np.array(errors_detrend))
				master_flags_detrend.append(np.array(flags_detrend))
			elif len(self.quarters) == 1:
				master_detrend = np.array(fluxes_detrend)
				master_error_detrend = np.array(errors_detrend)
				master_flags_detrend = np.array(flags_detrend)


		### this is the first initialization of the detrended fluxes.
		self.fluxes_detrend = master_detrend
		self.errors_detrend = master_error_detrend
		self.flags_detrend = master_flags_detrend 

		### before you overwrite the flags, compare them.

		if len(self.quarters) == 1:
			#final_flags = []
			#if type(self.flags_detrend) == list:
			#	self.flags_detrend = self.flags_detrend[0]
			#for sf, sfd in zip(self.flags, self.flags_detrend):
			#	final_flags.append(int(np.nanmax((sf,sfd))))
			final_flags = np.nanmax((self.flags, self.flags_detrend), axis=0)

		else:
			final_flags = []
			for qidx in np.arange(0,len(self.quarters),1):
				print('qidx = ', qidx)
				qfinal_flags = np.nanmax((self.flags[qidx], self.flags_detrend[qidx]), axis=0)
				print("qfinal_flags.shape = ", qfinal_flags.shape)
				final_flags.append(np.array(qfinal_flags))
			self.flags_detrend = final_flags

		if save_lc == 'y':
			lcfile = open(savepath+'/'+self.target+'_'+self.telescope+'_lightcurve.tsv', mode='w')
			if self.telescope.lower() == 'kepler' or self.telescope.lower() == 'k2':
				lcfile.write('BKJD\tfluxes\terrors\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
			elif self.telescope.lower() == 'tess':
				lcfile.write('BTJD\tfluxes\terrors\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
			elif self.telescope.lower() == 'user':
				lcfile.write('BJD\tfluxes\terrors\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
			### overwrite the existing file!
			self.fluxes_detrend = np.array(self.fluxes_detrend)
			self.errors_detrend = np.array(self.errors_detrend)
			self.flags_detrend = np.array(self.flags_detrend)
			lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = self.times, self.fluxes, self.errors, self.fluxes_detrend, self.errors_detrend, self.flags_detrend

			if len(self.quarters) > 1:
				for qidx in np.arange(0,len(self.quarters),1):
					print('qidx = ', qidx)
					qtq = self.quarters[qidx]
					qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags = lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx]

					for qt, qf, qe, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags):
						lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qfd)+'\t'+str(qed)+'\t'+str(qfl)+'\t'+str(qtq)+'\n')

					print('quarter written to file.')

			elif len(self.quarters) == 1:
				qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags = lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags
				
				if qtimes.shape[0] == 1:
					qtimes = qtimes[0]
				if qfluxes.shape[0] == 1:
					qfluxes = qfluxes[0]
				if qerrors.shape == 1:
					qerrors = qerrors[0]

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

				#print("qtimes = ", qtimes)
				#print('qfluxes = ', qfluxes)
				#print('qerrors = ', qerrors)
				#print('qfluxes_detrend = ', qfluxes_detrend)
				#print('qerrors_detrend = ', qerrors_detrend)
				#print('qflags = ', qflags)

				for qt, qf, qe, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags):
					lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qfd)+'\t'+str(qed)+'\t'+str(qfl)+'\t'+str(self.quarters[0])+'\n')

			lcfile.close()

		if self.telescope.lower() != 'user':
			self.get_neighbors(save_to_file='y')

		### finally, run get_neighbors() so that the neighbors will be appended to the end of the file.



	### FITTING!

	def fit(self, custom_param_dict=None, fitter='multinest', modelcode='LUNA', skip_ntqs='y', model='M', nlive=1000, nwalkers=100, nsteps=10000, resume=True, folded=False):
		### optional values for code are "multinest" and "emcee"
		#if type(params) != dict:
		#	raise Exception("'params' must be a dictionary, with strings as the keys and priors for the values.")

		### FOUR MODELS MAY BE RUN: a planet-only model (P) with no TTVs, a TTV model (T), freely fitting the transit times, 
		### a (Z) model, which gives the moon a mass but no radius, and an (M) model, which is a fully physical moon fit.

		if modelcode == "LUNA":
			folded = False ### you should not be fitting a moon model to a folded light curve!

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
			self.param_uber_dict['tau0'] = ['uniform', (self.fold_tau-0.1, self.fold_tau+0.1)]

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
			#nparamorig = 8  
			#nparam = 8
			nparamorig = 14 ### all these inputs must still be present, even if some of them are fixed at zero!
			nparam = 14
			nvars = 8  ### RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan
		elif model == 'T':
			#nparamorig = 8 ### RpRstar, rhostar, bplan, Pplan, tau0, q1, q2, rhoplan 
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
			param_labels.append(pkey)
			param_prior_forms.append(param_uber_dict[pkey][0])
			param_limit_tuple.append(param_uber_dict[pkey][1])

		self.param_labels = param_labels
		self.param_prior_forms = param_prior_forms 
		self.param_limit_tuple = param_limit_tuple


		for parlab, parprior, parlim in zip(param_labels, param_prior_forms, param_limit_tuple):
			print(parlab, parprior, parlim)


		if fitter == 'multinest':
			#mp_multinest(fit_times, fit_fluxes, fit_errors, param_labels=param_labels, param_prior_forms=param_prior_forms, param_limit_tuple=param_limit_tuple, nlive=nlive, targetID=self.target, modelcode=modelcode) ### outputs to a file
			mp_multinest(fit_times, fit_fluxes, fit_errors, param_dict=self.param_uber_dict, nlive=nlive, targetID=self.target, modelcode=modelcode, model=model, nparams=nvars)

		elif fitter == 'emcee':
			mp_emcee(fit_times, fit_fluxes, fit_errors, param_dict=self.param_uber_dict, nwalkers=nwalkers, nsteps=nsteps, targetID=self.target, modelcode=modelcode, model=model, resume=resume, nparams=nvars) ### outputs to a file

		### ONE FINAL STEP -- RESTORE DEFAULT VALUES (REMOVE tau0 = folded_tau) by initializing priors again.
		self.initialize_priors(modelcode=modelcode)		







	def genLS(self, show_plot = 'y', compute_fap='n', LSquarters=None):
		### this function generates a Lomb-Scargle Periodogram!
		LSperiods = []
		LSmaxpower_periods = []
		LSpowers = []
		LSfaps = []
		nquarters = len(self.quarters)

		if type(LSquarters) == type(None):
			LSquarters = self.quarters

		#for qidx in np.arange(0,nquarters,1):
		for qidx in np.arange(0,nquarters,1):
			this_quarter = self.quarters[qidx]
			if this_quarter not in LSquarters: ### use this in case you've only specified select quarters.
				continue

			print("processing LS for ", this_quarter)
			if nquarters != 1:
				qtimes, qfluxes, qerrors = self.times[qidx], self.fluxes[qidx], self.errors[qidx]
			else:
				qtimes, qfluxes, qerrors = self.times, self.fluxes, self.errors 
			maxperiod = 0.5 * (np.nanmax(qtimes) - np.nanmin(qtimes))
			minperiod = 0.5
			minfreq, maxfreq = 1/maxperiod, 1/minperiod
			qls = LombScargle(qtimes, qfluxes, qerrors)
			qfreq, qpower = qls.autopower(minimum_frequency=minfreq, maximum_frequency=maxfreq)
			qperiods = 1/qfreq
			if compute_fap == 'y':
				qfap = qls.false_alarm_probability(qpower.max(), method='bootstrap')
				probabilities = [0.1, 0.05, 0.01]
				quarter_FALs = qls.false_alarm_level(probabilities)

			if show_plot == 'y':
				#try:
				random_color = np.random.rand(3)
				plt.plot(qperiods[::-1], qpower[::-1], c=random_color)
				if compute_fap == 'y':
					plt.plot(qperiods[::-1], np.linspace(quarter_FALs[1], quarter_FALs[1], len(qperiods[::-1])), c=random_color)
				#except:
				#	plt.plot(qperiods[::-1], qpower[::-1])
				#	if compute_fap == 'y':
				#		plt.plot(qperiods[::-1], np.linspace(quarter_FALs[1], quarter_FALs[1], len(qperiods[::-1])))

			LSperiods.append(qperiods)
			max_power_period = qperiods[np.nanargmax(qpower)]
			LSmaxpower_periods.append(max_power_period)
			LSpowers.append(qpower)
			if compute_fap == 'y':
				LSfaps.append(qfap)

		if show_plot == 'y':
			plt.xscale('log')
			#plt.xlim(np.nanmin(qperiods), np.nanmax(qperiods))
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



	def find_transit_quarters(self, locate_neighbor='n'):
		self.get_properties(locate_neighbor=locate_neighbor)
		quarter_transit_dict = {}
		print('self.quarters = ', self.quarters)
		if len(self.quarters) > 1:
			print('option 1.')
			for qidx,quarter in enumerate(self.quarters):
				print('quarter = ', quarter)
				skip_this_quarter = 'n'
				quarter_times = self.times[qidx]
				### it's ok to do np.array(np.array([1])), doesn't nest them.
				### the code below throws an error when quarter_times is a float, so convert it to an array!
				desired_type = type(np.array([1]))
				if (type(quarter_times) == desired_type) == False:
					skip_this_quarter = 'y'
					continue

				if skip_this_quarter == 'y':
					continue

				if len(quarter_times) < 10: ### 5 hours.
					print('quarter has less than 10 data points.')
					continue

				#except:
				#	print('something goofy with this quarter. Skipping.')
				#	traceback.print_exc()
				#	skip_this_quarter = 'y'

				#if skip_this_quarter == 'y':
				#	continue

				quarter_times = np.array(quarter_times)

				quarter_transit = 'n'
				for tau in self.taus:
					if (tau >= np.nanmin(quarter_times)) and (tau <= np.nanmax(quarter_times)):
						### there's a transit in this quarter!
						quarter_transit = 'y'
						break

		elif len(self.quarters) == 1:
			print('option 2.')
			quarter_times = self.times
			quarter_transit = 'n'
			for tau in self.taus:
				if (tau >=- np.nanmin(quarter_times)) and (tau <= np.nanmax(quarter_times)):
					quarter_transit = 'y'
					break 

			if quarter_transit == 'n':
				quarter_transit_dict[0] = 'n'
			else:
				quarter_transit_dict[0] = 'y'
				
		self.quarter_transit_dict = quarter_transit_dict 


	def find_aliases(self):
		target_aliases = []
		alias_search_results = Simbad.query_objectids(self.target)
		for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
			target_alias = alias_search_results[alidx][0]
			target_aliases.append(target_alias)
		self.aliases = np.array(target_aliases)




	def examine_TPF(self, quarters=None, time_lims=None, detrend='y', mask_idxs=None):
		if type(quarters) == type(None):
			quarters = self.quarters 
		tpf_examiner(self.target, quarters=quarters, Tdur=self.duration_days, time_lims=time_lims, detrend=detrend, mask_idxs=mask_idxs, find_alias='y')





	def get_coords(self):
		targetID = self.targetID
		try:
			print(self.RA, self.Dec)
		except:
			if (self.telescope.lower() == 'user'):
				target_name = targetID
				target_in_simbad = 'n'

			elif (self.telescope.lower() == 'kepler'):
				if self.target_type == 'koi':
					target_name = "KOI-"+str(targetID)
				elif self.target_type == 'planet':
					target_name = 'Kepler-'+str(targetID)
				elif self.target_type == 'kic':
					target_name = "KIC "+str(targetID)

				try:
					simbad_query = Simbad.query_object(self.target)[0]
					self.RA = simbad_query['RA']
					self.Dec = simbad_query['DEC']
					target_in_simbad = 'y'
				except:
					target_in_simbad = 'n'
					try:
						self.RA = mast_data['ra'][self.mast_rowidx]
						self.Dec = mast_data['dec'][self.mast_rowidx]
					except:
						pass

			elif (self.telescope.lower() == 'k2'):
				target_name = "K2-"+str(targetID)
				try:
					simbad_query = Simbad.query_object(self.target)[0]
					self.RA = simbad_query['RA']
					self.Dec = simbad_query['DEC']
					target_in_simbad = 'y'
				except:
					target_in_simbad = 'n'
					try:
						self.RA = np.array(exofop_data['RA'][self.exofop_rowidx])[0]
						self.Dec = np.array(exofop_data["DEC"][self.exofop_rowidx])[0]
					except:
						pass


			elif (self.telescope.lower() == 'tess'):
				try:
					simbad_query = Simbad.query_object(self.target)[0] ### this will look for any object!
					target_name = self.target
					self.RA = simbad_query['RA']
					self.Dec = simbad_query['DEC']
					target_in_simbad = 'y'

				except:
					if self.target_type.lower() == 'toi':
						try:
							simbad_query = Simbad.query_object('TOI-'+self.target)[0]
							target_name = 'TOI-'+str(self.target) ### only update if this worked!
							self.RA = simbad_query['RA']
							self.Dec = simbad_query['DEC']
							target_in_simbad = 'y'
						except:
							target_in_simbad = 'n'
							try:
								self.RA = np.array(exofop_data['RA'][self.exofop_rowidx])[0]
								self.Dec = np.array(exofop_data['Dec'][self.exofop_rowidx])[0]
							except:
								pass

					elif self.target_type.lower() == 'tic':
						try:
							simbad_query = Simbad.query_object('TIC '+self.target)[0] 
							target_name = 'TIC '+str(self.target) ### only update if this worked!
							self.RA = simbad_query['RA']
							self.Dec = simbad_query['DEC']
							#target_in_simbad = 'y'
						except:
							#target_in_simbad = 'n'
							try:
								self.RA = np.array(exofop_data['RA'][self.exofop_rowidx])[0]
								self.Dec = np.array(exofop_data['Dec'][self.exofop_rowidx])[0]
							except:
								pass
				


	def fold(self, detrended='y'):
		### this method will phase fold your light curve. 
		### first tau in the time series:
		try:
			#first_tau = self.tau0
			first_tau = self.taus[0]
		except:
			self.get_properties()
			first_tau = self.taus[0]

		ftidx = 0
		while first_tau < np.nanmin(np.hstack(self.times)):
			ftidx += 1
			first_tau = self.taus[ftidx]

		fold_times = ((((np.hstack(self.times) - first_tau - 0.5*self.period) % self.period) / self.period)) ### yields the remainder!
		fold_times = fold_times - 0.5

		if detrended == 'y':
			fold_fluxes = np.hstack(self.fluxes_detrend)
			fold_errors = np.hstack(self.errors_detrend)
		else:
			fold_fluxes = np.hstack(self.fluxes)
			fold_errors = np.hstack(self.errors)

		self.fold_times = fold_times
		self.fold_fluxes = fold_fluxes
		self.fold_errors = fold_errors
		#self.fold_tau = fold_first_tau





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



	def find_planet_row(self, row_known='n'):
		download_new = 'n'
		if row_known == 'y':
			mast_rowidx = self.mast_rowidx
			if self.telescope.lower() != "kepler":
				exofop_rowidx = self.exofop_rowidx
		try:
			if self.telescope.lower() == 'user':
				download_new = 'n'

			elif (self.telescope.lower() == 'kepler'):
				filecreated_time = os.path.getctime(moonpydir+'/cumkois.txt')
				#filecreated_time2 = os.path.getctime(moonpydir+'/cfop_targets.csv')

			elif (self.telescope.lower() == 'k2'):
				filecreated_time1 = os.path.getctime(moonpydir+'/cumk2ois.txt')
				filecreated_time2 = os.path.getctime(moonpydir+'/exofop_targets.csv')
				filecreated_time = np.nanmin((filecreated_time1, filecreated_time2))

			elif (self.telescope.lower() == 'tess'):
				filecreated_time1 = os.path.getctime(moonpydir+'/exofop_toilists.pipe')
				filecreated_time2 = os.path.getctime(moonpydir+'/confirmed_planets.txt')
				filecreated_time = np.nanmin((filecreated_time1, filecreated_time2))

			current_time = time.time()
			if (current_time - filecreated_time) > 86400: ### the file is more than a day old.
				download_new = 'y'
		except:
			download_new = 'y'


		if download_new == 'y':
			### download a new version!
			try:
				if (self.telescope.lower() == 'kepler'):
					os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,ra,dec&order=kepoi_name&format=ascii" -O "'+moonpydir+'/cumkois.txt"')
					os.system('wget --tries=1 '+moonpydir+'/cfop_targets.csv "https://exofop.ipac.caltech.edu/kepler/download_summary_csv.php?sort=koi"')

					### TRY THIS IN THE COMMAND LINE.
					#wget --tries=5 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_prad,koi_prad_err1,koi_prad_err2,ra,dec&order=dec&format=ascii" -O "cumkois.txt"

				elif (self.telescope.lower() == 'k2'):
					os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=k2candidates&select=epic_name,epic_candname,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandep,pl_trandeperr1,pl_trandeperr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_eqt,pl_eqterr1,pl_eqterr2,pl_rade,pl_radeerr1,pl_radeerr2,st_rad,st_raderr1,st_raderr2,ra,dec&order=epic_name&format=ascii" -O "'+moonpydir+'/cumk2ois.txt"')
					os.system('wget --tries=1 '+moonpydir+'/exofop_targets.csv "https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All&sort=target"')

					### TRY THIS IN THE COMMAND LINE:
					#wget --tries=5 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=k2candidates&select=epic_name,epic_candname,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandep,pl_trandeperr1,pl_trandeperr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_eqt,pl_eqterr1,pl_eqterr2,pl_rade,pl_radeerr1,pl_radeerr2,st_rad,st_raderr1,st_raderr2,ra,dec&order=dec&format=ascii" -O "cumk2ois.txt"
					
					#wget --tries=5 exofop_targets.csv "https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All&sort=target"

				elif (self.telescope.lower() == 'tess'):
					os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_hostname,pl_letter,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_rade,pl_radeerr1,pl_radeerr2,ra,dec&order=pl_hostname&format=ascii" -O "'+moonpydir+'/confirmed_planets.txt"')
					os.system('wget --tries=1 '+moonpydir+'/exofop_toilists.pipe "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"')
			
					### TRY THIS ON THE COMMAND LINE:
					#wget --tries=5 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_hostname,pl_letter,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_rade,pl_radeerr1,pl_radeerr2,ra,dec&order=dec&format=ascii" -O "confirmed_planets.txt"
					
					#wget --tries=5 exofop_toilists.pipe "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"
			
			except:
				print("ATTEMPT TO DOWNLOAD THE LATEST PLANET CATALOGS FAILED (maybe you're not online?).")

		### USER INPUT HANDLING
		if self.telescope.lower() == 'user':
			mast_data = ascii.read('cumkois.txt')
			mast_columns = mast_data.columns
			exofop_data = pandas.read_csv('cfop_targets.csv', header=18)
			exofop_columns = exofop_data.columns
			row_known = 'y'
			self.mast_rowidx = np.nan 
			self.exofop_rowidx = np.nan 

		### KEPLER HANDLING 
		elif (self.telescope.lower() == 'kepler'):
			mast_data = ascii.read('cumkois.txt')
			mast_columns = mast_data.columns
			try:
				exofop_data = pandas.read_csv('cfop_targets.csv', header=18)
				exofop_columns = exofop_data.columns
			except:
				pass

			if row_known == 'n':
				if str(self.target).lower().startswith('kepler'):
					target_letter = str(self.target)[-1]
					if ' ' in self.target: ### already in the correct format, with a space between the letter.
						NEA_targetname = self.target
					else: #### there isn't a space, so add it!
						NEA_targetname = self.target[:-1]+' '+target_letter
					mast_rowidx = np.where(mast_data['kepler_name'] == NEA_targetname)[0]


				elif str(self.target).lower().startswith('kic'):
					NEA_targetname = int(self.target[4:])
					mast_rowidx = np.where(mast_data['kepid'] == NEA_targetname)[0]

				elif str(self.target).lower().startswith('koi'):
					NEA_targetname = str(self.target[4:])
					if len(NEA_targetname) == 7: ### of the form 5084.01
						NEA_targetname = 'K0'+str(NEA_targetname)
					elif len(NEA_targetname) == 6: ### of the form 163.01
						NEA_targetname = 'K00'+str(NEA_targetname)
					elif len(NEA_targetname) == 5: ### of the form 23.01
						NEA_targetname = 'K000'+str(NEA_targetname)
					elif len(NEA_targetname) == 4: ### of the form 1.01
						NEA_targetname = 'K0000'+str(NEA_targetname)
					mast_rowidx = np.where(mast_data['kepoi_name'] == NEA_targetname)[0]
					self.NEA_targetname = NEA_targetname


				else: ### was observed by Kepler but you don't know it's KOI/KIC or Kepler name:
					mast_data = ascii.read('confirmed_planets.txt') ### overwrite before! won't be found in cumkois!
					try:
						float(self.target[-1]) ### if this works, query the mast_data['pl_hostname'] because you end with a number!
						NEA_taretname = self.target 
						mast_rowidx = np.where(mast_data['pl_hostname'] == self.target)[0]

					except:
						### implies the last thing is a number
						target_letter = self.target[-1]
						if ' ' in self.target:
							NEA_targetname = self.target
						else:
							NEA_targetname = self.target[:-1]+' '+target_letter
						mast_rowidx = np.where(mast_data['pl_name'] == NEA_targetname)[0]

				self.mast_rowidx = mast_rowidx 


		### K2 HANDLING 
		elif (self.telescope.lower() == 'k2'):
			mast_data = ascii.read('cumk2ois.txt')
			mast_columns = mast_data.columns
			exofop_data = pandas.read_csv('exofop_targets.csv', header=10)
			exofop_columns = exofop_data.columns

			if row_known == 'n':
				if str(self.target).startswith('K2-'):
					target_letter = str(self.target[-1])
					if ' ' in self.target:
						NEA_targetname = self.target
					else:
						NEA_targetname = str(self.target[:-1])+' '+str(self.target[-1])
					mast_rowidx = np.where(mast_data['pl_name'] == NEA_targetname)[0]
					exofop_rowidx = np.nan

					print("number of rows matching this description = ", len(mast_rowidx))

				elif str(self.target).startswith('EPIC'):
					if self.target.startswith('EPIC'):
						NEA_targetname = self.target[4:] ### just the numbers!
						if (NEA_targetname.startswith(' ')) or (NEA_targetname.startswith('-')):
							NEA_targetname = NEA_targetname[1:]
					try:
						exofop_rowidx = np.where(exofop_data['EPIC ID'] == NEA_targetname)[0]
					except:
						print('unable to extract the exofop_rowidx')
						exofop_rowidx = np.nan 

					try:
						mast_rowidx = np.where(mast_data['epic_name'] == NEA_targetname)[0]
					except:
						print('unable to extract the mast_rowidx')
						mast_rowidx = np.nan


		### TESS HANDLING 
		elif (self.telescope.lower() == 'tess'):
			mast_data = ascii.read('confirmed_planets.txt')
			mast_columns = mast_data.columns
			exofop_data = pandas.read_csv('exofop_toilists.pipe', delimiter='|')
			exofop_columns = exofop_data.columns

			if row_known == 'n':
				try:
					float(self.target[-1]) 
					### if the except statement isn't triggered, the last value is a number! Therefore
					NEA_targetname = str(self.target)
				except: 
					target_letter = str(self.target[-1])
					### implies the last value is a letter.
					if ' ' in self.target:
						NEA_targetname = str(self.target) ### we want a space in here.
					else: ### put one in!
						NEA_targetname = str(self.target[:-1])+' '+str(target_letter)

				### FOR TESS TARGETS THE EXOFOP TABLE IS LIKELY BETTER THAN THE NEA table.
				if NEA_targetname.lower().startswith('tic'):
					### search the "TIC ID" column
					ticnum = NEA_targetname[3:]
					if (ticnum.startswith(' ')) or (ticnum.startswith('-')):
						ticnum = ticnum[1:]
					exofop_rowidx = np.where(exofop_data['TIC ID'] == int(ticnum))[0]
					print('looking for '+str(ticnum)+' in the exofop database.')
					#print('rowidx (pre-clean) = ', rowidx)
					if len(exofop_rowidx) > 1:
						exofop_rowidx = rowidx[0]
					print("exofop_rowidx = ", exofop_rowidx)
					mast_rowidx = np.nan 

				elif NEA_targetname.lower().startswith("toi"):
					### search the TOI column 
					toinum = NEA_targetname[3:]
					if (toinum.startswith(' ')) or (toinum.startswith('-')):
						toinum = toinum[1:]
					print('looking for '+str(toinum)+' in the exofop database.')
					exofop_rowidx = np.where(exofop_data['TOI'] == float(toinum))[0]
					print('exofop_rowidx = ', exofop_rowidx)
					mast_rowidx = np.nan 

				else:
					### in this case, you have to look in the mast list!
					mast_rowidx = np.where(mast_data['pl_hostname'] == NEA_targetname)[0] ### in "confirmed" planets.
					nfound = 0
					exofop_rowidxs = []
					for idx in np.arange(0,len(exofop_data['Comments']),1): ### have to use a loop here.
						try:
							float(self.target[-1]) ### if this works, last digit is a number, so use self.target
							if (self.target in exofop_data['Comments'][idx]):
								exofop_rowidx = idx
								exofop_rowidxs.append(exofop_rowidx)
								print('found a planet in CFOP!, exofop_rowidx, TIC = ', exofop_rowidx, exofop_data['TIC ID'][idx])
								#print('TIC = ', exofop_data['TIC ID'][exofop_rowidx])
								nfound += 1
						except:
							### last digit is a letter, so get rid of it when searching!
							try:
								if (self.target[:-1] in np.array(exofop_data['Comments'])[idx]):
									exofop_rowidx = idx
									exofop_rowidxs.append(exofop_rowidx)
									print('found a planet in CFOP!, exofop_rowidx, TIC = ', exofop_rowidx, exofop_data['TIC ID'][idx])
									#print("TIC = ", exofop_data['TIC ID'][exofop_rowidx])
									nfound += 1
							except:
								pass

					print('# CFOP entries: ', nfound)

					if nfound > 1:
						print('found multiple possible entries for '+self.target)
						for exfidx in exofop_rowidxs:
							print("idx, TIC, comment = ", exfidx, exofop_data['TIC ID'][exfidx], exofop_data['Comments'][exfidx])
						print(' ')
						exofop_rowidx = int(input('Which is the correct index? ')) 

					elif nfound == 1:
						exofop_rowidx = exofop_rowidxs[0]

					elif nfound == 0:
						### try aliases!
						exofop_rowidx = np.nan
						try:
							for alias in self.aliases:
								mast_rowidx = np.where(mast_data['pl_hostname'] == alias)[0]
								if len(mast_rowidx) != 0:
									break
						except:
							pass

		if row_known == 'n':
			try:
				self.exofop_rowidx = exofop_rowidx
			except:
				pass

			try:
				self.mast_rowidx = mast_rowidx
			except:
				pass

		#print('rowidx = ', rowidx)
		if self.telescope.lower() == 'user':
			try:
				return mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname
			except:
				return mast_rowidx, exofop_rowidx, mast_data, exofop_data, self.target 			

		elif self.telescope.lower() == 'k2' or self.telescope.lower() == 'tess':
			try:
				return mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname
			except:
				return mast_rowidx, exofop_rowidx, mast_data, exofop_data, self.target 
		elif self.telescope.lower() == 'kepler':
			try:
				return mast_rowidx, mast_data, NEA_targetname
			except:
				return mast_rowidx, mast_data, self.target







	def get_properties(self, locate_neighbor='n'):
		print("calling 'get_properties()...")
		if self.telescope.lower() == 'user':
			target_period, target_period_uperr, target_period_lowerr = self.period, 0, 0
			target_tau0, target_tau0_uperr, target_tau0_lowerr = self.tau0, 0, 0
			target_impact, target_impact_uperr, target_impact_lowerr = self.impact, 0, 0
			target_duration, target_duration_uperr, target_duration_lowerr = self.duration_hours, 0, 0
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = self.rprstar, 0, 0
			target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = self.sma_AU, 0, 0
			target_insol, target_insol_uperr, target_insol_lowerr = np.nan, 0, 0
			target_rp, target_rp_uperr, target_rp_lowerr = self.rp_rearth, 0, 0

		elif self.telescope.lower() == 'k2' or self.telescope.lower() == 'tess':
			if self.newlc == 'y':
				mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row(row_known='y')
			elif self.newlc == 'n':
				mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row(row_known='n')

			print('mast_rowidx = ', mast_rowidx)
			print('exofop_rowidx = ', exofop_rowidx)
		else:
			if self.newlc == 'y':
				mast_rowidx, mast_data, NEA_targetname = self.find_planet_row(row_known='y')
			elif self.newlc == 'n': 
				mast_rowidx, mast_data, NEA_targetname = self.find_planet_row(row_known='n')	
			print('mast_rowidx = ', mast_rowidx)

		### now with the rowidx we can access the other properties we want!
		if (self.telescope.lower() == 'kepler'):
			target_period, target_period_uperr, target_period_lowerr = mast_data['koi_period'][mast_rowidx], mast_data['koi_period_err1'][mast_rowidx], mast_data['koi_period_err2'][mast_rowidx]
			target_tau0, target_tau0_uperr, target_tau0_lowerr = mast_data['koi_time0bk'][mast_rowidx], mast_data['koi_time0bk_err1'][mast_rowidx], mast_data['koi_time0bk_err2'][mast_rowidx]
			target_impact, target_impact_uperr, target_impact_lowerr = mast_data['koi_impact'][mast_rowidx], mast_data['koi_impact_err1'][mast_rowidx], mast_data['koi_impact_err2'][mast_rowidx]
			target_duration, target_duration_uperr, target_duration_lowerr = mast_data['koi_duration'][mast_rowidx], mast_data['koi_duration_err1'][mast_rowidx], mast_data['koi_duration_err2'][mast_rowidx]
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = mast_data['koi_ror'][mast_rowidx], mast_data['koi_ror_err1'][mast_rowidx], mast_data['koi_ror_err2'][mast_rowidx]
			target_rp, target_rp_uperr, target_rp_lowerr = mast_data['koi_prad'][mast_rowidx], mast_data['koi_prad_err1'][mast_rowidx], mast_data['koi_prad_err2'][mast_rowidx]
			target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = mast_data['koi_sma'][mast_rowidx], mast_data['koi_sma_err1'][mast_rowidx], mast_data['koi_sma_err2'][mast_rowidx]
			target_insol, target_insol_uperr, target_insol_lowerr = mast_data['koi_insol'][mast_rowidx], mast_data['koi_insol_err1'][mast_rowidx], mast_data['koi_insol_err2'][mast_rowidx]
			target_ldm1, target_ldm2 = mast_data['koi_ldm_coeff1'][mast_rowidx], mast_data['koi_ldm_coeff2'][mast_rowidx]
			target_eccen, target_eccen_uperr, target_eccen_lowerr = mast_data['koi_eccen'][mast_rowidx], mast_data['koi_eccen_err1'][mast_rowidx], mast_data['koi_eccen_err2'][mast_rowidx]
			target_longp, target_longp_uperr, target_longp_lowerr = mast_data['koi_longp'][mast_rowidx], mast_data['koi_longp_err1'][mast_rowidx], mast_data['koi_longp_err2'][mast_rowidx]
			target_incl, target_incl_uperr, target_incl_lowerr = mast_data['koi_incl'][mast_rowidx], mast_data['koi_incl_err1'][mast_rowidx], mast_data['koi_incl_err2'][mast_rowidx]	



		elif (self.telescope.lower() == 'k2'):
			target_period, target_period_uperr, target_period_lowerr = np.nanmedian(mast_data['pl_orbper'][mast_rowidx]), np.nanmedian(mast_data['pl_orbpererr1'][mast_rowidx]), np.nanmedian(mast_data['pl_orbpererr2'][mast_rowidx])
			target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nanmedian(mast_data['pl_tranmid'][mast_rowidx]), np.nanmedian(mast_data['pl_tranmiderr1'][mast_rowidx]), np.nanmedian(mast_data['pl_tranmiderr2'][mast_rowidx])
			target_impact, target_impact_uperr, target_impact_lowerr = np.nanmedian(mast_data['pl_imppar'][mast_rowidx]), np.nanmedian(mast_data['pl_impparerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_impparerr2'][mast_rowidx])
			target_duration, target_duration_uperr, target_duration_lowerr = np.nanmedian(mast_data['pl_trandur'][mast_rowidx]), np.nanmedian(mast_data['pl_trandurerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_trandurerr2'][mast_rowidx])
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.nanmedian(mast_data['pl_ratror'][mast_rowidx]), np.nanmedian(mast_data['pl_ratrorerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_ratrorerr2'][mast_rowidx])
			target_rp, target_rp_uperr, target_rp_lowerr = np.nanmedian(mast_data['pl_rade'][mast_rowidx]), np.nanmedian(mast_data['pl_radeerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_radeerr2'][mast_rowidx])
			target_a_rstar, target_a_rstar_uperr, target_a_rstar_lowerr = np.nanmedian(mast_data['pl_ratdor'][mast_rowidx]), np.nanmedian(mast_data['pl_ratdorerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_ratdorerr2'][mast_rowidx])
			target_Teq, target_Teq_uperr, target_Teq_lowerr = np.nanmedian(mast_data['pl_eqt'][mast_rowidx]), np.nanmedian(mast_data['pl_eqterr1'][mast_rowidx]), np.nanmedian(mast_data['pl_eqterr2'][mast_rowidx])


		elif (self.telescope.lower() == "tess"):
			#if (self.target.startswith('TIC')) or (self.target.startswith("TOI")):
			if (np.isfinite(self.exofop_rowidx) == True) or (np.isfinite(self.exofop_rowidx) == np.array([True])):
				print('searching for target parameters in the exofop database.')
				print(' ')
				target_period, target_period_uperr, target_period_lowerr = np.array(exofop_data['Period (days)'])[exofop_rowidx], np.array(exofop_data['Epoch (BJD) err'])[exofop_rowidx], np.array(exofop_data['Epoch (BJD) err'])[exofop_rowidx]
				target_impact, target_impact_uperr, target_impact_lowerr = np.nan, np.nan, np.nan 
				target_duration, target_duration_uperr, target_duration_lowerr = np.array(exofop_data['Duration (hours)'])[exofop_rowidx], np.array(exofop_data['Duration (hours) err'])[exofop_rowidx], np.array(exofop_data['Duration (hours) err'])[exofop_rowidx]
				target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.sqrt(1e-6*np.array(exofop_data['Depth (ppm)'])[exofop_rowidx]), np.sqrt(1e-6*np.array(exofop_data['Depth (ppm) err'])[exofop_rowidx]), np.sqrt(1e-6*np.array(exofop_data['Depth (ppm) err'])[exofop_rowidx])
				target_rp, target_rp_uperr, target_rp_lowerr = np.array(exofop_data['Planet Radius (R_Earth)'])[exofop_rowidx], np.array(exofop_data['Planet Radius (R_Earth) err'])[exofop_rowidx], np.array(exofop_data['Planet Radius (R_Earth) err'])[exofop_rowidx]
				target_tau0 = np.nan

			elif np.isfinite(self.mast_rowidx) == True:
				print('search for target parameters in the mast database.')
				print(' ')
				target_period, target_period_uperr, target_period_lowerr = np.nanmedian(mast_data['pl_orbper'][mast_rowidx]), np.nanmedian(mast_data['pl_orbpererr1'][mast_rowidx]), np.nanmedian(mast_data['pl_orbpererr2'][mast_rowidx])
				target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nanmedian(mast_data['pl_tranmid'][mast_rowidx]), np.nanmedian(mast_data['pl_tranmiderr1'][mast_rowidx]), np.nanmedian(mast_data['pl_tranmiderr2'][mast_rowidx])
				target_impact, target_impact_uperr, target_impact_lowerr = np.nanmedian(mast_data['pl_imppar'][mast_rowidx]), np.nanmedian(mast_data['pl_impparerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_impparerr2'][mast_rowidx])
				target_duration, target_duration_uperr, target_duration_lowerr = np.nanmedian(mast_data['pl_trandur'][mast_rowidx]), np.nanmedian(mast_data['pl_trandurerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_trandurerr2'][mast_rowidx])
				target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.nanmedian(mast_data['pl_ratror'][mast_rowidx]), np.nanmedian(mast_data['pl_ratrorerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_ratrorerr2'][mast_rowidx])
				target_rp, target_rp_uperr, target_rp_lowerr = np.nanmedian(mast_data['pl_rade'][mast_rowidx]), np.nanmedian(mast_data['pl_radeerr1'][mast_rowidx]), np.nanmedian(mast_data['pl_radeerr2'][mast_rowidx])
			else:
				pass

		### update properties!
		try:
			self.period = float(target_period)
		except:
			print('target_period = ', target_period)
			raise Exception('a problem was encountered converting target_period to a float.')
		try:
			self.period_err = (float(target_period_lowerr), float(target_period_uperr))
		except:
			self.period_err = (np.nan, np.nan)
		if (self.telescope.lower() == 'tess') and (float(target_tau0) > 2454833):
			try:
				self.tau0 = float(target_tau0) - 2457000 ### native TESS offset value
			except:
				self.tau0 = np.nan
		else:
			self.tau0 = float(target_tau0)
		try:
			self.tau0_err = (float(target_tau0_lowerr), float(target_tau0_uperr))
		except:
			self.tau0_err = (np.nan, np.nan)

		if np.isfinite(target_impact):
			self.impact = float(target_impact)
			self.impact_err = (float(target_impact_lowerr), float(target_impact_uperr))
		else:
			print('impact parameter is NaN. Setting = 0.')
			self.impact = 0
			self.impact_err = (0,0)
		self.rprstar = float(target_rprstar)
		self.rprstar_err = (float(target_rprstar_lowerr), float(target_rprstar_uperr))
		self.rp_rearth = float(target_rp) ### earth radii
		self.rp_rjup = float(target_rp) * (R_earth.value / R_jup.value)
		self.rp_meters = self.rp_rjup * eq_RJup
		self.rp_rearth_err = (float(target_rp_lowerr), float(target_rp_uperr))
		self.rstar_rsol = (float(target_rp) * (1/float(target_rprstar))) * (R_earth.value / R_sun.value)
		self.rstar_meters = self.rstar_rsol * eq_RSun
		self.depth = self.rprstar**2
		if self.telescope == 'kepler':
			self.ldm_a1 = float(target_ldm1)
			self.ldm_a2 = float(target_ldm2)
			self.q1, self.q2 = u1u2_to_q1q2(self.ldm_a1, self.ldm_a2)

		try:
			self.sma_AU = float(target_sma_AU) ### sma in AU
			self.sma_AU_err = (float(target_sma_AU_lowerr), float(target_sma_AU_uperr)) ### sma in AU
		except:
			try:
				target_sma_AU = (self.a_rstar * self.rstar_meters) / au.value 
				self.sma_AU = target_sma_AU 
			except:
				pass

		try:
			self.a_rstar = float(target_a_rstar) ### sma/Rstar
			self.a_rstar_err = (float(target_a_rstar_lowerr), float(target_a_rstar_uperr)) ### sma/Rstar 
		except:
			### try to convert sma_AU into A/Rstar
			try:
				target_a_rstar = (self.sma_AU * au.value) / self.rstar_meters
				self.a_rstar = target_a_rstar 
			except:
				pass

		try:
			self.insol = float(target_insol) ### insolution in units of Earth insolation
			self.insol_err = (float(target_insol_lowerr), float(target_insol_uperr))
		except:
			pass

		try:
			self.Teq = float(target_Teq) ### equilibrium temperature of the planet
			self.Teq_err = (float(target_Teq_lowerr), float(target_Teq_uperr))
		except:
			pass


		try:
			self.eccen = float(target_eccen) ### equilibrium temperature of the planet
			self.eccen_err = (float(target_eccen_lowerr), float(target_eccen_uperr))
		except:
			pass

		try:
			self.longp = float(target_longp) ### equilibrium temperature of the planet
			self.longp_err = (float(target_longp_lowerr), float(target_longp_uperr))
		except:
			pass


		try:
			self.incl = float(target_incl) ### equilibrium temperature of the planet
			self.incl_err = (float(target_incl_lowerr), float(target_incl_uperr))
		except:
			pass



		if np.isfinite(target_duration):
			self.duration_hours = float(target_duration)
			self.duration_hours_err = (float(target_duration_lowerr), float(target_duration_uperr))
			self.duration_days = float(target_duration)/24
			self.duration_days_err = (float(target_duration_lowerr)/24, float(target_duration_uperr)/24)
		else:
			print('transit duration is NaN. Calculating an estimate.')
			### calculate it!
			### estimate the mass of the star from its radius.
			target_mstar_msol_est = self.rstar_rsol ** 1.1364
			target_mstar_kg_est = MSun*target_mstar_msol_est
			target_est_sma = Kep3_afromp(self.period, target_mstar_kg_est, 1e-3*target_mstar_kg_est)
			print('target_est_sma = ', target_est_sma)
			target_est_Tdur = Tdur(self.period, self.rstar_meters, self.rp_meters, self.impact, target_est_sma)
			print('target_est_Tdur = ', target_est_Tdur)
			self.duration_days = target_est_Tdur
			self.duration_days_err = (-0.01*self.duration_days, 0.01*self.duration_days)
			self.duration_hours = self.duration_days * 24
			self.duration_hours_err = self.duration_days_err * 24 

		if locate_neighbor=='y':
			self.find_neighbors() ### new May 31st, 2019 -- identifies whether there are other known planets in the system!
		try:
			self.find_taus()
		except:
			traceback.print_exc()
			print("UNABLE TO CALCULATE transit times. You may not have downloaded any data.")
		###	identify in-transit times


	def find_taus(self):
		try:
			print("calling 'find_taus()'.")
			transit_midtimes = [self.tau0]
			while (transit_midtimes[-1] - self.period) > np.nanmin(np.hstack(self.times)):
				### the transit_midtime you just added isn't the first transit!
				transit_midtimes.append(transit_midtimes[-1] - self.period)
			transit_midtimes = np.sort(transit_midtimes).tolist()

			next_transit = transit_midtimes[-1]+self.period
			nquarters = len(self.quarters)
			if nquarters != 1:
				maxtime = np.nanmax(np.concatenate((self.times)))
			elif nquarters == 1:
				maxtime = np.nanmax(self.times)
			while next_transit < maxtime:
				transit_midtimes.append(next_transit)
				next_transit = transit_midtimes[-1]+self.period
			self.taus = np.array(transit_midtimes)
		except:
			traceback.print_exc()
			raise Exception('an exception was raised while calling find_taus().')


	def mystery_solver(self, tau0, period, duration_hours, neighbor_tau0=None, neighbor_period=None, neighbor_duration_hours=None, neighbor_name='None'):
		### this function will generate essential 
		self.tau0 = tau0
		self.period = period
		self.duration_hours = duration_hours
		self.duration_days = duration_hours / 24
		self.find_taus()
		### calculate transit times for this target... assuming linear ephemeris!
		transit_times = []
		nquarters = len(self.quarters)
		for qidx in np.arange(0,nquarters,1):
			if nquarters != 1:
				qtimes = self.times[qidx]
			elif nquarters == 1:
				qtimes = self.times
			for qtime in qtimes:
				for tau in self.taus:
					if (qtime >= tau - 0.5*self.duration_days) and (qtime <= tau + 0.5*self.duration_days):
						transit_times.append(qtime)
		transit_times = np.unique(np.array(transit_times))
		self.all_transit_times = transit_times 
		self.neighbors = []






	def find_neighbors(self, is_neighbor='n'):
		print('is_neighbor = ', is_neighbor)
		if is_neighbor == 'y':
			row_known = 'n'
		else:
			row_known = 'y'

		if self.telescope.lower() == 'user':
			pass 

		elif self.telescope.lower() == 'k2' or self.telescope.lower() == 'tess':
			mast_rowidx, exofop_rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row(row_known=row_known)
			#check_mast_rows = np.arange(0,len(mast_data['pl_name']),1)
			#check_exofop_rows = np.arange(0,len(exofop_data['TOI']),1)
			if mast_rowidx < 10:
				check_mast_rows = np.arange(0,mast_rowidx+11,1)
			else:
				check_mast_rows = np.arange(mast_rowidx-10,mast_rowidx+10,1)

			if exofop_rowidx < 10:
				check_exofop_rows = np.arange(0,exofop_rowidx+11,1)
			else:
				check_exofop_rows = np.arange(exofop_rowidx-10,exofop_rowidx+11,1)
			print('mast_rowidx, exofop_rowidx = ', mast_rowidx, exofop_rowidx)
		
		else:
			mast_rowidx, mast_data, NEA_targetname = self.find_planet_row(row_known=row_known)
			#check_mast_rows = np.arange(0,len(mast_data['kepoi_name']),1)
			if mast_rowidx < 10:
				check_mast_rows = np.arange(0,mast_rowidx+11,1)
			else:
				check_mast_rows = np.arange(mast_rowidx-10,mast_rowidx+10,1)

			print('mast_rowidx = ', mast_rowidx)



		neighbor_rows = []
		neighbor_targets = []

		if self.target.lower().startswith('usr'):
			pass

		elif self.target.lower().startswith('koi'):
			print("looking for neighbors in MAST for this KOI.")

			for cr in check_mast_rows:
				if cr <= len(mast_data['kepoi_name']) - 1:
					if (np.array(mast_data['kepoi_name'])[cr][:-1] == self.NEA_targetname[:-1]) and (cr != mast_rowidx):
						print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['kepoi_name'][cr]))
						neighbor = str(mast_data['kepoi_name'][cr])
						while neighbor.startswith('K') or neighbor.startswith('0'):
							neighbor = neighbor[1:]
						neighbor = 'KOI-'+str(neighbor)
						neighbor_rows.append(cr)
						neighbor_targets.append(neighbor)


		elif self.target.lower().startswith('kepler'):
			print("looking for neighbors in MAST for this Kepler target.")
			for cr in check_mast_rows:
				if cr <= len(mast_data['kepler_name']) - 1:
					if (np.array(mast_data['kepler_name'])[cr][7:-2] == self.target[7:-1]) and (cr != mast_rowidx):
						print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['kepler_name'][cr]))
						neighbor = str(mast_data['kepler_name'][cr])
						if ' ' in neighbor:
							neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
						neighbor_rows.append(cr)
						neighbor_targets.append(neighbor)

		elif self.target.lower().startswith('k2'):
			print("looking for neighbors in MAST for this K2 target.")
			for cr in check_mast_rows:
				if cr <= len(mast_data['pl_name']) - 1:
					if (np.array(mast_data['pl_name'])[cr][3:-2] == NEA_targetname[3:-2]) and (np.array(mast_data['pl_name'])[cr] != NEA_targetname):
						print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['pl_name'][cr]))
						neighbor = str(mast_data['pl_name'][cr])
						if ' ' in neighbor:
							neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
						neighbor_rows.append(cr)
						neighbor_targets.append(neighbor)		

		else:
			if (self.telescope.lower() == "tess"):
				if self.target.lower().startswith("tic"):
					print("looking for neighbors in exofop for this TIC target.")
					for cr in check_exofop_rows:
						if cr <= len(exofop_data['TIC ID']) - 1:
							if (str(np.array(exofop_data['TIC ID'])[cr])[3:] == str(NEA_targetname)[3:]) and (cr != exofop_rowidx):
								print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(exofop_data['TIC ID'][cr]))
								neighbor = str(exofop_data['TIC ID'][cr])
								if ' ' in neighbor:
									neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
								neighbor_rows.append(cr)
								neighbor_targets.append('TIC '+str(neighbor))

				elif self.target.lower().startswith('toi'):
					print('looking for neighbors in exofop for this TOI.')
					for cr in check_exofop_rows:
						if cr <= len(exofop_data['TOI']) - 1:
							if (str(np.array(exofop_data['TOI'])[cr])[:-2] == str(NEA_targetname)[4:-2]) and (cr != exofop_rowidx):
								print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(exofop_data['TOI'][cr]))
								neighbor = str(exofop_data['TOI'][cr])
								if ' ' in neighbor:
									neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
								neighbor_rows.append(cr)
								neighbor_targets.append('TOI-'+str(neighbor))						

				else:
					try:
						print('looking fore neighbors in MAST for this TESS target.')
						for cr in check_mast_rows:
							if cr <= len(mast_data['pl_name']) - 1:
								if (np.array(mast_data['pl_name'])[cr][3:-2] == NEA_targetname[3:-2]) and (np.array(mast_data['pl_name'])[cr] != NEA_targetname):
									print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['pl_name'][cr]))
									neighbor = str(mast_data['pl_name'][cr])
									if ' ' in neighbor:
										neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
									neighbor_rows.append(cr)
									neighbor_targets.append(neighbor)		
					except:
						pass							

		self.neighbors = neighbor_targets



	def get_neighbors(self, clobber_lc='y', save_to_file='y'):
		### this function will download grab all the information about your target's neighbors.
		if self.telescope.lower() == 'user':
			neighbor_dict = {}

		try:
			print(self.neighbor_dict.keys())
			neighbor_dict = self.neighbor_dict
			### indicates you're dealing with the target. Therefore:
			is_neighbor='n'
		except:
			is_neighbor='y'
			### you need to generate the neighbor_dict
			neighbor_dict = {}
			for neighbor in self.neighbors:
				###
				try:
					if (self.telescope.lower() == 'kepler'):
						if neighbor.lower().startswith('kepler'):
							neighbor_key = 'k'+str(neighbor[7:])
						else:
							### for now this is just a KOI!
							neighbor_key = neighbor
							while neighbor_key.startswith('K') or neighbor_key.startswith('0'):
								neighbor_key = neighbor_key[1:]
							neighbor_key = 'k'+str(neighbor_key)
					elif (self.telescope.lower() == 'k2'):
						neighbor_key = 'k2_'+str(neighbor_key[3:])

					elif (self.telescope.lower() == 'tess'):
						neighbor_key = neighbor 

					### now generate the LC object -- note that this will download duplicate light curves!
					### for now, just delete them after downloading... should come up with a better solution.
					print('calling MoonpyLC for neighbor: ', neighbor)

					neighbor_dict[neighbor_key] = MoonpyLC(targetID=neighbor, is_neighbor=is_neighbor, clobber='n')
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

			for nt in neighbor_taus:
				ntidxs = np.where((np.hstack(self.times) >= (nt - 0.5*neighbor_dur)) & (np.hstack(self.times) <= (nt + 0.5*neighbor_dur)))[0]
				for ntidx in ntidxs:
					neighbor_transit_idxs.append(ntidx)
					neighbor_transit_IDs.append(neighbor)

		### include the target here!
		target_taus = self.taus
		target_dur = self.duration_days
		for tt in target_taus:
			ttidxs = np.where((np.hstack(self.times) >= (tt - 0.5*target_dur)) & (np.hstack(self.times) <= (tt + 0.5*target_dur)))[0]
			for ttidx in ttidxs:
				neighbor_transit_idxs.append(ttidx)
				neighbor_transit_IDs.append(self.target)

		#neighbor_transit_idxs = np.hstack(neighbor_transit_idxs)
		#neighbor_transit_IDs = np.hstack(neighbor_transit_IDs)

		neighbor_transit_idxs = np.array(neighbor_transit_idxs)
		neighbor_transit_IDs = np.array(neighbor_transit_IDs)

		for ntidx in np.unique(neighbor_transit_idxs):
			### find all transiting planets that have this index
			all_neighbors = neighbor_transit_IDs[np.where(neighbor_transit_idxs == ntidx)[0]]
			final_neighbor_IDs.append(list(all_neighbors))


		neighbor_transit_idxs = np.unique(neighbor_transit_idxs)

		neighbor_transit_times = []
		neighbor_transit_list = []
		neighbor_transit_ID = []

		if save_to_file == 'y':
			try:
				lcfile = open(savepath+'/'+self.target+"_"+self.telescope+"_lightcurve.tsv", mode='r')
			except:
				lcfile = open(savepath+'/'+self.target+'_lightcurve.tsv', mode='r') ### for older files.
			lcfile_new = open(savepath+"/"+self.target+"_"+self.telescope+"_lc_temp.tsv", mode='w')

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
			#os.system('mv '+savepath+'/'+self.target+'_lc_temp.tsv '+savepath+'/'+self.target+'_lightcurve.tsv')
			os.system('mv '+savepath+'/'+self.target+'_'+self.telescope+'_lc_temp.tsv '+savepath+'/'+self.target+'_'+self.telescope+'_lightcurve.tsv')

		self.all_transit_times = np.array(neighbor_transit_times)
		self.all_transit_list = neighbor_transit_list
		self.all_transit_IDs = neighbor_transit_ID 



	def prep_for_CNN(self, save_lc='y', window=6, cnn_len=493, exclude_neighbors='y', flag_neighbors='y', show_plot='n', extra_path_info=None, cnnlc_path=moonpydir+'/cnn_lcs'):
		### this function will produce an arrayy that's ready to be fed to a CNN for moon finding. 
		if os.path.exists(cnnlc_path) == False:
			os.system('mkdir '+cnnlc_path)

		localpath_list = []	

		for taunum,tau in enumerate(self.taus):
			cnn_min, cnn_max = tau-window, tau+window
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
							ntidxs = np.where((cnn_times >= nt-(0.5*self.neighbor_dict[neighbor].duration_days)) & (cnn_times <= nt+(0.5*self.neighbor_dict[neighbor].duration_days)))[0]
							flag_idxs.append(ntidxs)
				try:
					flag_idxs = np.unique(np.hstack(np.array(flag_idxs)))
				except:
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




	def find_TTVs(self, show_plot='n', yvar='OCmins', window=2):
		OC_mins = []
		OC_days = []
		OC_over_durs = []
		epochs = []

		nquarters = len(self.quarters)
		#for qidx,q in enumerate(self.quarters):
		for qidx in np.arange(0,quarters,1):
			if nquarters != 1:
				qtimes, qfluxes, qerrors = self.times[qidx], self.fluxes_detrend[qidx], self.errors_detrend[qidx]
			elif nquarters == 1:
				qtimes, qfluxes, qerrors = self.times, self.fluxes_detrend, self.errors_detrend 
			for epoch, tau in enumerate(self.taus):
				if tau >= np.nanmin(qtimes) and tau <= np.nanmax(qtimes):
					### this tau is in the quarter
					### grab times up to 3x the transit duration on either side of the the tau.
					transit_idxs = np.where((qtimes >= tau - (window*self.duration_days)) & (qtimes <= tau + (window*self.duration_days)))[0]
					transit_times = qtimes[transit_idxs]
					transit_fluxes = qfluxes[transit_idxs]
					numerator = np.nansum((1 - transit_fluxes) * transit_times)
					denominator = np.nansum(1 - transit_fluxes)
					tmid = numerator / denominator
					OC_day = tmid - tau
					OC_over_dur = OC_day / self.duration_days
					OC_min = OC_day * (24 * 60)

					OC_mins.append(OC_min)
					OC_days.append(OC_day)
					OC_over_durs.append(OC_over_dur)
					epochs.append(epoch)
		OC_mins = np.array(OC_mins)
		OC_days = np.array(OC_days)
		OC_over_durs = np.array(OC_over_durs)
		epochs = np.array(epochs)
		OC_sig_min = np.nanstd(OC_mins)
		OC_sig_day = np.nanstd(OC_days)
		OC_sig_over_dur = OC_sig_day / self.duration_days

		self.OCs_min = OC_mins
		self.OCs_day = OC_days
		self.OCs_over_dur = OC_over_durs
		self.OC_sig_min = OC_sig_min
		self.OC_sig_day = OC_sig_day
		self.OC_sig_over_dur = OC_sig_over_dur

		if show_plot == 'y':
			if yvar == "OCmins":
				yvals = OC_mins
				ylab = 'O - C [minutes]'
			elif yvar == 'OCdays':
				yvals = OC_days
				ylab = 'O - C [days]'
			elif yvar == 'OCdurs':
				yvals = OC_over_durs
				ylab = '(O - C) / duration'

			plt.scatter(epochs, yvals, s=10, facecolor='LightCoral', edgecolor='k')
			plt.xlabel('Epoch')
			plt.ylabel(ylab)
			plt.show()









	###############################
	### VISUALIZATION FUNCTIONS ###
	###############################


	def plot_lc(self, facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', folded='n', include_flagged='n', detrended='y', show_errors='n', show_neighbors='n', time_format='native', pltshow='y'):
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
			if time_format == 'native':
				plot_stitched_times = stitched_times
			elif time_format == 'bjd':
				if self.telescope == 'kepler':
					plot_stitched_times = stitched_times + 2454833
				elif self.telescope == 'tess':
					plot_stitched_times = stitched_times + 2457000

			plt.scatter(plot_stitched_times, stitched_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			if show_errors == 'y':
				plt.errorbar(plot_stitched_times, stitched_fluxes, yerr=stitched_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')

			if show_neighbors == 'y':
				### this will highlight all the other transits for the neighbors (if any)
				neighbors = self.neighbor_dict.keys()

				for neighbor in neighbors:
					neighbor_taus = self.neighbor_dict[neighbor].taus 
					neighbor_dur = self.neighbor_dict[neighbor].duration_days 

					neighbor_transit_idxs = []
					for nt in neighbor_taus:
						ntidxs = np.where((stitched_times >= (nt - 0.5*neighbor_dur)) & (stitched_times <= (nt + 0.5*neighbor_dur)))[0]
						neighbor_transit_idxs.append(ntidxs)
					neighbor_transit_idxs = np.hstack(neighbor_transit_idxs)
					#plt.scatter(stitched_times[neighbor_transit_idxs], stitched_fluxes[neighbor_transit_idxs], facecolors='g', s=10, marker='x')
					plt.scatter(plot_stitched_times[neighbor_transit_idxs], stitched_fluxes[neighbor_transit_idxs], s=10, marker='x', label=neighbor)
			

				### PLOT THE TARGET TRANSITS TOO!
				target_taus = self.taus 
				target_dur = self.duration_days
				target_transit_idxs = []
				for tt in target_taus:
					ttidxs = np.where((stitched_times >= (tt - 0.5*target_dur)) & (stitched_times <= (tt + 0.5*target_dur)))[0]
					target_transit_idxs.append(ttidxs)
				target_transit_idxs = np.hstack(target_transit_idxs)
				plt.scatter(plot_stitched_times[target_transit_idxs], stitched_fluxes[target_transit_idxs], s=10, marker='x', color='Indigo', label='target')

			plt.legend()					


		elif folded == 'y':
			try:
				self.fold(detrended=detrended)
			except:
				self.get_properties(locate_neighbor='n')
				self.fold()
			plt.scatter(self.fold_times, self.fold_fluxes, facecolors=facecolor, edgecolors=edgecolor, s=10, zorder=1)
			if show_errors == 'y':
				plt.errorbar(self.fold_times, self.fold_fluxes, yerr=self.fold_errors, ecolor='k', zorder=0, alpha=0.5, fmt='none')
		#plt.xlabel('BKJD')
		if (self.telescope.lower() == 'kepler') or (self.telescope.lower() == 'k2'):
			if folded=='y':
				plt.xlabel('Phase')
			else:
				plt.xlabel('BKJD')
		elif (self.telescope.lower() == 'tess'):
			if folded=='y':
				plt.xlabel('Phase')
			else:
				plt.xlabel('BTJD')
		plt.ylabel('Flux')
		try:
			plt.title(str(self.target))
		except:
			pass

		if pltshow == 'y':	
			plt.show()
		else:
			pass



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


