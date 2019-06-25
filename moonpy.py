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

	def __init__(self, targetID=None, target_type=None, quarters='all', telescope=None, RA=None, Dec=None, coord_format='degrees', search_radius=5, lc_format='pdc', remove_flagged='y', sc=False, ffi='n', save_lc='y', load_lc='n', is_neighbor='n', clobber=None):

		if (load_lc == 'n') and (clobber == None):
			### check to see if a file already exists!
			if os.path.exists(savepath+'/'+str(targetID)+'_lightcurve.tsv'):
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


		if load_lc == 'y':
			if self.target.startswith('K2') or self.target.startswith('k2'):
				self.telescope = "k2"
			elif self.target.startswith('Kepler') or self.target.startswith("kepler") or self.target.startswith('KIC') or self.target.startswith('kic') or self.target.startswith('KOI') or self.target.startswith('koi'):
				self.telescope = "kepler"
			elif self.target.startswith('TIC') or self.target.startswith('TOI'):
				self.telescope = 'tess'
			else:
				telescope = input('Please specify the telescope: ')
				self.telescope = telescope 

			try:
				pandafile = pandas.read_csv(savepath+'/'+target_name+'_lightcurve.tsv', delimiter='\t')
				ptimes = np.array(pandafile['BKJD'])
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



		### HANDLING FOR DOWNLOADING A FRESH LIGHT CURVE.
		if (load_lc=='n') and (targetID != None) and (telescope != None):
			### implies you've selected a target you want to download.
			if (telescope == 'kepler') or (telescope=="Kepler") or (telescope == 'KEPLER'):
				rowidx, mast_data, NEA_targetname = self.find_planet_row()
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, type=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)

			elif (telescope == 'K2') or (telescope == 'k2'):
				rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row()
				lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = kplr_target_download(targetID, type=target_type, quarters=quarters, telescope=telescope, lc_format=lc_format, sc=sc)

			elif (telescope == 'tess') or (telescope == "Tess") or (telescope == "TESS"):
				rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row()

				if ffi == 'y':
					lc_times, lc_fluxes, lc_errors = eleanor_target_download(targetID, lc_format=lc_format, sectors=quarters)

				elif ffi == 'n':
					if self.target.startswith("TOI"):
						ticnum = int(np.array(exofop_data['TIC ID'])[rowidx])
						ticname = 'TIC '+str(ticnum)
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(ticname)	
						print('lc_times.shape = ', lc_times.shape)
					else:
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_target_download(targetID)
						print('lc_times.shape = ', lc_times.shape)

			self.telescope = telescope

			if (telescope == 'kepler') or (telescope == 'Kepler') or (telescope == "KEPLER"):
				if target_type == 'koi':
					target_name = "KOI-"+str(targetID)
				elif target_type == 'planet':
					target_name = 'Kepler-'+str(targetID)
				elif target_type == 'kic':
					target_name = "KIC "+str(targetID)

				try:
					simbad_query = Simbad.query_object(self.target)[0]
					self.RA = simbad_query['RA']
					self.Dec = simbad_query['DEC']
					target_in_simbad = 'y'
				except:
					target_in_simbad = 'n'
					try:
						self.RA = mast_data['ra'][rowidx]
						self.Dec = mast_data['dec'][rowidx]
					except:
						pass


			elif (telescope == 'k2') or (telescope == "K2"):
				target_name = "K2-"+str(targetID)
				try:
					simbad_query = Simbad.query_object(self.target)[0]
					self.RA = simbad_query['RA']
					self.Dec = simbad_query['DEC']
					target_in_simbad = 'y'
				except:
					target_in_simbad = 'n'
					try:
						self.RA = np.array(exofop_data['RA'][rowidx])[0]
						self.Dec = np.array(exofop_data["DEC"][rowidx])[0]
					except:
						pass


			elif (telescope == 'TESS') or (telescope == 'tess') or (telescope == 'Tess'):
				try:
					simbad_query = Simbad.query_object(self.target)[0]
					target_name = self.target
					self.RA = simbad_query['RA']
					self.Dec = simbad_query['DEC']
					target_in_simbad = 'y'

				except:
					if target_type == 'toi':
						try:
							simbad_query = Simbad.query_object('TOI-'+self.target)[0]
							target_name = 'TOI-'+str(self.target) ### only update if this worked!
							self.RA = simbad_query['RA']
							self.Dec = simbad_query['DEC']
							target_in_simbad = 'y'
						except:
							target_in_simbad = 'n'
							try:
								self.RA = np.array(exofop_data['RA'][rowidx])[0]
								self.Dec = np.array(exofop_data['Dec'][rowidx])[0]
							except:
								pass

					elif target_type == 'tic':
						try:
							simbad_query = Simbad.query_object('TIC '+self.target)[0] 
							target_name = 'TIC '+str(self.target) ### only update if this worked!
							self.RA = simbad_query['RA']
							self.Dec = simbad_query['DEC']
							#target_in_simbad = 'y'
						except:
							#target_in_simbad = 'n'
							try:
								self.RA = np.array(exofop_data['RA'][rowidx])[0]
								self.Dec = np.array(exofop_data['Dec'][rowidx])[0]
							except:
								pass

				try:
					self.find_aliases()
				except:
					print('ALIASES COULD NOT BE FOUND.')


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
					try:
						lc_times, lc_fluxes, lc_errors, lc_flags, lc_quarters = tess_coord_download(RA, Dec)
					except:
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
			if len(lc_quarters) > 1:
				timesort = np.argsort(lc_times[qidx])

				if lc_fluxes_detrend.shape == 0:
					lc_fluxes_detrend = lc_fluxes
				if lc_errors_detrend.shape == 0:
					lc_errors_detrend = lc_errors 
				lc_times[qidx], lc_fluxes[qidx], lc_errors[qidx], lc_fluxes_detrend[qidx], lc_errors_detrend[qidx], lc_flags[qidx] = lc_times[qidx][timesort], lc_fluxes[qidx][timesort], lc_errors[qidx][timesort], lc_fluxes_detrend[qidx][timesort], lc_errors_detrend[qidx][timesort], lc_flags[qidx][timesort]
			
			elif len(lc_quarters) == 1:
				timesort = np.argsort(lc_times)
				lc_times, lc_fluxes, lc_errors, lc_fluxes_detrend, lc_errors_detrend, lc_flags = lc_times[timesort], lc_fluxes[timesort], lc_errors[timesort], lc_fluxes_detrend[timesort], lc_errors_detrend[timesort], lc_flags[timesort]


		self.times = lc_times
		self.fluxes = lc_fluxes
		self.errors = lc_errors
		self.flags = lc_flags
		self.quarters = lc_quarters

		### calculate which quarters have transits!
		try:
			#self.find_transit_quarters()
			if is_neighbor == 'y':
				#self.get_properties(locate_neighbor='n')
				self.find_transit_quarters(locate_neighbor='n') ### get_properties() is called first thing!
			elif is_neighbor == 'n':
				self.find_transit_quarters(locate_neighbor='y') ### ditto above.
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


		if save_lc == 'y':
			### write to a file!
			lcfile = open(savepath+'/'+str(target_name)+'_lightcurve.tsv', mode='w')
			lcfile.write('BKJD\tfluxes\terrors\tflags\tquarter\n')
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







	### DETRENDING!

	def detrend(self, dmeth='cofiam', save_lc='y', mask_transits='y', mask_neighbors='y', skip_ntqs='n', kernel=None, max_degree=30):
		exceptions_raised = 'n'

		if mask_neighbors == 'y':
			mask_transits = 'y' 
		### optional values for method are "cofiam", "untrendy", "medfilt"
		### EACH QUARTER SHOULD BE DETRENDED INDIVIDUALLY!

		try:
			self.duration_hours ### tests whether you've called get_properties yet.
		except:
			self.get_properties(locate_neighbor='n')

		master_detrend, master_error_detrend, master_flags_detrend = [], [], []

		nquarters = len(self.quarters)
		for qidx in np.arange(0,nquarters,1):
			skip_quarter = 'n'
			print('quarter = ', self.quarters[qidx])
			if nquarters != 1:
				dtimes, dfluxes, derrors = self.times[qidx], self.fluxes[qidx], self.errors[qidx]
			elif nquarters == 1:
				dtimes, dfluxes, derrors = self.times, self.fluxes, self.errors
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
				try:
					quarter_transit_taus = self.taus[((self.taus > np.nanmin(dtimes)) & (self.taus < np.nanmax(dtimes)))]
				except:
					self.find_transit_quarters()
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
					neighbor_transit_idxs = []
					for ntt in self.all_transit_times:
						if ntt >= np.nanmin(dtimes) and ntt <= np.nanmax(dtimes):
							neighbor_transit_idxs.append(np.where(ntt == dtimes)[0])
					mask_transit_idxs.append(neighbor_transit_idxs)


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
						flags_detrend = np.linspace(0,0,len(fluxes_detrend))

					elif dmeth == 'untrendy':
						print("UNTRENDY ISN'T REALLY SUPPORTED RIGHT NOW, SORRY!")
						fluxes_detrend, errors_detrend = untrendy_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs)
						flags_detrend = np.linspace(0,0,len(fluxes_detrend))

					elif dmeth == 'george':
						print("GEORGE ISN'T REALLY SUPPORTED RIGHT NOW, SORRY!")
						fluxes_detrend, errors_detrend = george_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, mask_idxs=mask_transit_idxs)
						flags_detrend = np.linspace(0,0,len(fluxes_detrend))

					elif dmeth == 'medfilt':
						print("MEDIAN FILTERING HAS YET TO BE TESTED EXTENSIVELY. BEWARE!")
						fluxes_detrend, errors_detrend = medfilt_detrend(dtimes, dfluxes, derrors, telescope=self.telescope, size=kernel, mask_idxs=mask_transit_idxs)
						flags_detrend = np.linspace(0,0,len(fluxes_detrend))

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


			master_detrend.append(np.array(fluxes_detrend))
			master_error_detrend.append(np.array(errors_detrend))
			master_flags_detrend.append(np.array(flags_detrend))


		### this is the first initialization of the detrended fluxes.
		self.fluxes_detrend = master_detrend
		self.errors_detrend = master_error_detrend
		self.flags_detrend = master_flags_detrend 

		### before you overwrite the flags, compare them.
		final_flags = []
		for qidx in np.arange(0,len(self.quarters),1):
			qfinal_flags = []

			if type(self.flags[0]) == np.int32:
				self.flags = [self.flags]
			for sf,sfd in zip(self.flags[qidx], self.flags_detrend[qidx]):
				qfinal_flags.append(int(np.nanmax((sf,sfd))))


			final_flags.append(np.array(qfinal_flags))
		self.flags_detrend = final_flags

		if save_lc == 'y':
			lcfile = open(savepath+'/'+self.target+'_lightcurve.tsv', mode='w')
			lcfile.write('BKJD\tfluxes\terrors\tfluxes_detrended\terrors_detrended\tflags\tquarter\n')
			### overwrite the existing file!
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

				try:
					if qflags.shape == 1:
						qflags = qflags[0]
				except:
					if type(qflags) == list:
						qflags = qflags[0]


				for qt, qf, qe, qfd, qed, qfl in zip(qtimes, qfluxes, qerrors, qfluxes_detrend, qerrors_detrend, qflags):
					lcfile.write(str(qt)+'\t'+str(qf)+'\t'+str(qe)+'\t'+str(qfd)+'\t'+str(qed)+'\t'+str(qfl)+'\t'+str(self.quarters[0])+'\n')

			lcfile.close()

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

		self.get_properties(locate_neighbor='n')
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







	def genLS(self, show_plot = 'y'):
		### this function generates a Lomb-Scargle Periodogram!
		LSperiods = []
		LSpowers = []
		LSfaps = []
		nquarters = len(self.quarters)

		for qidx in np.arange(0,nquarters,1):
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


	def find_transit_quarters(self, locate_neighbor='n'):
		self.get_properties(locate_neighbor=locate_neighbor)
		quarter_transit_dict = {}
		if len(self.quarters > 1):
			for qidx,quarter in enumerate(self.quarters):
				quarter_times = self.times[qidx]
				quarter_transit = 'n'
				for tau in self.taus:
					if (tau >= np.nanmin(quarter_times)) and (tau <= np.nanmax(quarter_times)):
						### there's a transit in this quarter!
						quarter_transit = 'y'
						break

		elif len(self.quarters == 1):
			quarter_times = self.times
			quarter_transit = 'n'
			for tau in self.taus:
				if (tau >=- np.nanmin(quarter_times)) and (tau <= np.nanmax(quarter_times)):
					quarter_transit = 'y'
					break 

			if quarter_transit == 'n':
				quarter_transit_dict[quarter] = 'n'
			else:
				quarter_transit_dict[quarter] = 'y'
				
		self.quarter_transit_dict = quarter_transit_dict 


	def find_aliases(self):
		target_aliases = []
		alias_search_results = Simbad.query_objectids(self.target)
		for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
			target_alias = alias_search_results[alidx][0]
			target_aliases.append(target_alias)
		self.aliases = np.array(target_aliases)



	def fold(self, detrended='y'):
		### this method will phase fold your light curve. 
		### first tau in the time series:
		self.get_properties(locate_neighbor='n')
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



	def find_planet_row(self):
		download_new = 'n'
		try:
			if (self.telescope == 'kepler') or (self.telescope == "Kepler"):
				filecreated_time1 = os.path.getctime(moonpydir+'/cumkois.txt')
				filecreated_time2 = os.path.getctime(moonpydir+'/cfop_targets.csv')

			elif (self.telescope == 'k2') or (self.telescope == 'K2'):
				filecreated_time1 = os.path.getctime(moonpydir+'/cumk2ois.txt')
				filecreated_time2 = os.path.getctime(moonpydir+'/exofop_targets.csv')

			elif (self.telescope == 'tess') or (self.telescope == 'TESS') or (self.telescope == "Tess"):
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
			if (self.telescope == 'kepler') or (self.telescope == "Kepler"):
				os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_prad,koi_prad_err1,koi_prad_err2,ra,dec&order=dec&format=ascii" -O "'+moonpydir+'/cumkois.txt"')
				os.system('wget -O '+moonpydir+'/cfop_targets.csv "https://exofop.ipac.caltech.edu/kepler/download_summary_csv.php?sort=koi"')

			elif (self.telescope == 'k2') or (self.telescope == "K2"):
				os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=k2candidates&select=epic_name,epic_candname,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandep,pl_trandeperr1,pl_trandeperr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_eqt,pl_eqterr1,pl_eqterr2,pl_rade,pl_radeerr1,pl_radeerr2,st_rad,st_raderr1,st_raderr2,ra,dec&order=dec&format=ascii" -O "'+moonpydir+'/cumk2ois.txt"')
				os.system('wget -O '+moonpydir+'/exofop_targets.csv "https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All&sort=target"')

			elif (self.telescope == 'TESS') or (self.telescope == 'Tess') or (self.telescope == 'tess'):
				os.system('wget "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_hostname,pl_letter,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_rade,pl_radeerr1,pl_radeerr2,ra,dec&order=dec&format=ascii" -O "'+moonpydir+'/confirmed_planets.txt"')
				os.system('wget -O '+moonpydir+'/exofop_toilists.pipe "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=pipe"')


		### find by KICID, KOI number of planet!
		if (self.telescope == 'kepler') or (self.telescope == "Kepler"):
			mast_data = ascii.read('cumkois.txt')
			mast_columns = mast_data.columns
			try:
				exofop_data = pandas.read_csv('cfop_targets.csv', header=18)
				exofop_columns = exofop_data.columns
			except:
				pass

		elif (self.telescope == 'k2') or (self.telescope == 'K2'):
			mast_data = ascii.read('cumk2ois.txt')
			mast_columns = mast_data.columns
			exofop_data = pandas.read_csv('exofop_targets.csv', header=10)
			exofop_columns = exofop_data.columns

		elif (self.telescope == 'TESS') or (self.telescope == 'Tess') or (self.telescope == 'tess'):
			mast_data = ascii.read('confirmed_planets.txt')
			mast_columns = mast_data.columns
			exofop_data = pandas.read_csv('exofop_toilists.pipe', delimiter='|')
			exofop_columns = exofop_data.columns



		### KEPLER PRIMARY MISSION HANDLING
		if (self.telescope == 'Kepler') or (self.telescope == 'kepler'):
			if str(self.target).startswith('Kepler-'):
				target_letter = str(self.target)[-1]
				if ' ' in self.target: ### already in the correct format, with a space between the letter.
					NEA_targetname = self.target
				else: #### there isn't a space, so add it!
					NEA_targetname = self.target[:-1]+' '+target_letter
				rowidx = np.where(mast_data['kepler_name'] == NEA_targetname)[0]

			elif str(self.target).startswith('KIC'):
				NEA_targetname = int(self.target[4:])
				rowidx = np.where(mast_data['kepid'] == NEA_targetname)[0]

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
				rowidx = np.where(mast_data['kepoi_name'] == NEA_targetname)[0]



		### K2 MISSION HANDLING
		elif (self.telescope == 'k2') or (self.telescope == 'K2'):
			if str(self.target).startswith('K2-'):
				target_letter = str(self.target[-1])
				if ' ' in self.target:
					NEA_targetname = self.target
				else:
					NEA_targetname = str(self.target[:-1])+' '+str(self.target[-1])
				rowidx = np.where(mast_data['pl_name'] == NEA_targetname)[0]

				print("number of rows matching this description = ", len(rowidx))

			elif str(self.target).startswith('EPIC'):
				if self.target.startswith('EPIC'):
					NEA_targetname = self.target[4:] ### just the numbers!
					if (NEA_targetname.startswith(' ')) or (NEA_targetname.startswith('-')):
						NEA_targetname = NEA_targetname[1:]
				try:
					rowidx = np.where(exofop_data['EPIC ID'] == NEA_targetname)[0]
				except:
					print('exofop_data.columns = ')
					print(exofop_data.columns)	
					time.sleep(60)

		### TESS MISSION HANDLING
		elif (self.telescope == 'TESS') or (self.telescope == "Tess") or (self.telescope == 'tess'):
			try:
				float(self.target[-1])
				### if the except statement isn't triggered, the last value is a number! Therefore
				NEA_targetname = str(self.target)
			except:
				target_letter = str(self.target[-1])
				### implies the last value is a letter.
				if ' ' in self.target:
					NEA_targetname = str(self.target)
				else:
					NEA_targetname = str(self.target[:-1])+' '+str(target_letter)

			### FOR TESS TARGETS THE EXOFOP TABLE IS LIKELY BETTER THAN THE NEA table.
			if NEA_targetname.startswith('TIC'):
				### search the "TIC ID" column
				ticnum = NEA_targetname[3:]
				if (ticnum.startswith(' ')) or (ticnum.startswith('-')):
					ticnum = ticnum[1:]
				rowidx = np.where(exofop_data['TIC ID'] == int(ticnum))[0]
				print('looking for '+str(ticnum)+' in the exofop database.')
				#print('rowidx (pre-clean) = ', rowidx)
				if len(rowidx) > 1:
					rowidx = rowidx[0]
				print("rowidx = ", rowidx)


			elif NEA_targetname.startswith("TOI"):
				toinum = NEA_targetname[3:]
				if (toinum.startswith(' ')) or (toinum.startswith('-')):
					toinum = toinum[1:]
				print('looking for '+str(toinum)+' in the exofop database.')
				rowidx = np.where(exofop_data['TOI'] == float(toinum))[0]
				print('rowidx = ', rowidx)



			else:
				rowidx = np.where(mast_data['pl_name'] == NEA_targetname)[0]
				if len(rowidx) == 0:
					rowidx = np.where(mast_data['pl_hostname'] == NEA_targetname)[0]
					### try aliases!
					if len(rowidx) == 0: ### if you STILL haven't found it.
						try:
							for alias in self.aliases:
								rowidx = np.where(mast_data['pl_hostname'] == alias)[0]
								if len(rowidx) != 0:
									break
						except:
							pass



		#print('rowidx = ', rowidx)
		if self.telescope == 'K2' or self.telescope == 'k2' or self.telescope == 'TESS' or self.telescope == 'Tess' or self.telescope == 'tess':
			return rowidx, mast_data, exofop_data, NEA_targetname
		else:
			return rowidx, mast_data, NEA_targetname




	def get_properties(self, locate_neighbor='n'):
		if self.telescope == 'K2' or self.telescope == 'k2' or self.telescope == 'TESS' or self.telescope == 'Tess' or self.telescope == 'tess':
			rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row()
		else:
			rowidx, mast_data, NEA_targetname = self.find_planet_row()

		### now with the rowidx we can access the other properties we want!
		if (self.telescope == 'Kepler') or (self.telescope == 'kepler'):
			target_period, target_period_uperr, target_period_lowerr = mast_data['koi_period'][rowidx], mast_data['koi_period_err1'][rowidx], mast_data['koi_period_err2'][rowidx]
			target_tau0, target_tau0_uperr, target_tau0_lowerr = mast_data['koi_time0bk'][rowidx], mast_data['koi_time0bk_err1'][rowidx], mast_data['koi_time0bk_err2'][rowidx]
			target_impact, target_impact_uperr, target_impact_lowerr = mast_data['koi_impact'][rowidx], mast_data['koi_impact_err1'][rowidx], mast_data['koi_impact_err2'][rowidx]
			target_duration, target_duration_uperr, target_duration_lowerr = mast_data['koi_duration'][rowidx], mast_data['koi_duration_err1'][rowidx], mast_data['koi_duration_err2'][rowidx]
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = mast_data['koi_ror'][rowidx], mast_data['koi_ror_err1'][rowidx], mast_data['koi_ror_err2'][rowidx]
			target_rp, target_rp_uperr, target_rp_lowerr = mast_data['koi_prad'][rowidx], mast_data['koi_prad_err1'][rowidx], mast_data['koi_prad_err2'][rowidx]
			target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = mast_data['koi_sma'][rowidx], mast_data['koi_sma_err1'][rowidx], mast_data['koi_sma_err2'][rowidx]
			target_insol, target_insol_uperr, target_insol_lowerr = mast_data['koi_insol'][rowidx], mast_data['koi_insol_err1'][rowidx], mast_data['koi_insol_err2'][rowidx]


		elif (self.telescope == 'k2') or (self.telescope == "K2"):
			target_period, target_period_uperr, target_period_lowerr = np.nanmedian(mast_data['pl_orbper'][rowidx]), np.nanmedian(mast_data['pl_orbpererr1'][rowidx]), np.nanmedian(mast_data['pl_orbpererr2'][rowidx])
			target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nanmedian(mast_data['pl_tranmid'][rowidx]), np.nanmedian(mast_data['pl_tranmiderr1'][rowidx]), np.nanmedian(mast_data['pl_tranmiderr2'][rowidx])
			target_impact, target_impact_uperr, target_impact_lowerr = np.nanmedian(mast_data['pl_imppar'][rowidx]), np.nanmedian(mast_data['pl_impparerr1'][rowidx]), np.nanmedian(mast_data['pl_impparerr2'][rowidx])
			target_duration, target_duration_uperr, target_duration_lowerr = np.nanmedian(mast_data['pl_trandur'][rowidx]), np.nanmedian(mast_data['pl_trandurerr1'][rowidx]), np.nanmedian(mast_data['pl_trandurerr2'][rowidx])
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.nanmedian(mast_data['pl_ratror'][rowidx]), np.nanmedian(mast_data['pl_ratrorerr1'][rowidx]), np.nanmedian(mast_data['pl_ratrorerr2'][rowidx])
			target_rp, target_rp_uperr, target_rp_lowerr = np.nanmedian(mast_data['pl_rade'][rowidx]), np.nanmedian(mast_data['pl_radeerr1'][rowidx]), np.nanmedian(mast_data['pl_radeerr2'][rowidx])
			target_a_rstar, target_a_rstar_uperr, target_a_rstar_lowerr = np.nanmedian(mast_data['pl_ratdor'][rowidx]), np.nanmedian(mast_data['pl_ratdorerr1'][rowidx]), np.nanmedian(mast_data['pl_ratdorerr2'][rowidx])
			target_Teq, target_Teq_uperr, target_Teq_lowerr = np.nanmedian(mast_data['pl_eqt'][rowidx]), np.nanmedian(mast_data['pl_eqterr1'][rowidx]), np.nanmedian(mast_data['pl_eqterr2'][rowidx])


		elif (self.telescope == "TESS") or (self.telescope == "Tess") or (self.telescope == 'tess'):
			if (self.target.startswith('TIC')) or (self.target.startswith("TOI")):
				print('searching for target parameters in the exofop database.')
				target_period, target_period_uperr, target_period_lowerr = np.array(exofop_data['Period (days)'])[rowidx], np.array(exofop_data['Period (days) err'])[rowidx], np.array(exofop_data['Period (days) err'])[rowidx]
				target_tau0, target_tau0_uperr, target_tau0_lowerr = np.array(exofop_data['Epoch (BJD)'])[rowidx], np.array(exofop_data['Epoch (BJD) err'])[rowidx], np.array(exofop_data['Epoch (BJD) err'])[rowidx]
				target_impact, target_impact_uperr, target_impact_lowerr = np.nan, np.nan, np.nan 
				target_duration, target_duration_uperr, target_duration_lowerr = np.array(exofop_data['Duration (hours)'])[rowidx], np.array(exofop_data['Duration (hours) err'])[rowidx], np.array(exofop_data['Duration (hours) err'])[rowidx]
				target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.sqrt(1e-6*np.array(exofop_data['Depth (ppm)'])[rowidx]), np.sqrt(1e-6*np.array(exofop_data['Depth (ppm) err'])[rowidx]), np.sqrt(1e-6*np.array(exofop_data['Depth (ppm) err'])[rowidx])
				target_rp, target_rp_uperr, target_rp_lowerr = np.array(exofop_data['Planet Radius (R_Earth)'])[rowidx], np.array(exofop_data['Planet Radius (R_Earth) err'])[rowidx], np.array(exofop_data['Planet Radius (R_Earth) err'])[rowidx]
			
			else:
				print('search for target parameters in the mast database.')
				target_period, target_period_uperr, target_period_lowerr = np.nanmedian(mast_data['pl_orbper'][rowidx]), np.nanmedian(mast_data['pl_orbpererr1'][rowidx]), np.nanmedian(mast_data['pl_orbpererr2'][rowidx])
				target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nanmedian(mast_data['pl_tranmid'][rowidx]), np.nanmedian(mast_data['pl_tranmiderr1'][rowidx]), np.nanmedian(mast_data['pl_tranmiderr2'][rowidx])
				target_impact, target_impact_uperr, target_impact_lowerr = np.nanmedian(mast_data['pl_imppar'][rowidx]), np.nanmedian(mast_data['pl_impparerr1'][rowidx]), np.nanmedian(mast_data['pl_impparerr2'][rowidx])
				target_duration, target_duration_uperr, target_duration_lowerr = np.nanmedian(mast_data['pl_trandur'][rowidx]), np.nanmedian(mast_data['pl_trandurerr1'][rowidx]), np.nanmedian(mast_data['pl_trandurerr2'][rowidx])
				target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.nanmedian(mast_data['pl_ratror'][rowidx]), np.nanmedian(mast_data['pl_ratrorerr1'][rowidx]), np.nanmedian(mast_data['pl_ratrorerr2'][rowidx])
				target_rp, target_rp_uperr, target_rp_lowerr = np.nanmedian(mast_data['pl_rade'][rowidx]), np.nanmedian(mast_data['pl_radeerr1'][rowidx]), np.nanmedian(mast_data['pl_radeerr2'][rowidx])

		"""
		print('target_period = ', target_period)
		print('target_tau0 = ', target_tau0)
		print('target_impact = ', target_impact)
		print("target_duration = ", target_duration)
		print('target_rprstar = ', target_rprstar)
		print('target_rp = ', target_rp)
		"""


		### update properties!
		self.period = float(target_period)
		self.period_err = (float(target_period_lowerr), float(target_period_uperr))
		if (self.telescope == 'tess') and (float(target_tau0) > 2454833):
			self.tau0 = float(target_tau0) - 2457000 ### native TESS offset value
		else:
			self.tau0 = float(target_tau0)
		self.tau0_err = (float(target_tau0_lowerr), float(target_tau0_uperr))
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

		self.find_taus()
		###	identify in-transit times


	def find_taus(self):
		transit_midtimes = [self.tau0]
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
		print('calling find_neighbors. is_neighbor = ', is_neighbor)
		if self.telescope == 'K2' or self.telescope == 'k2' or self.telescope == 'TESS' or self.telescope == 'Tess' or self.telescope == 'tess':
			rowidx, mast_data, exofop_data, NEA_targetname = self.find_planet_row()
		else:
			rowidx, mast_data, NEA_targetname = self.find_planet_row()


		### NOW IDENTIFY WHETHER THERE ARE ANY PLANETS IN THE VICINITY OF THIS PLANET (rowidx+/-7) with similar names!
		if rowidx < 10:
			check_rows = np.arange(0,rowidx+11,1)
		else:
			check_rows = np.arange(rowidx-10,rowidx+11,1)
		neighbor_rows = []
		neighbor_targets = []

		if (self.target).startswith('KOI') or (self.target).startswith('koi'):
			for cr in check_rows:
				if (np.array(mast_data['kepoi_name'])[cr][:-1] == NEA_targetname[:-1]) and (cr != rowidx):
					print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['kepoi_name'][cr]))
					neighbor = str(mast_data['kepoi_name'][cr])
					while neighbor.startswith('K') or neighbor.startswith('0'):
						neighbor = neighbor[1:]
					neighbor = 'KOI-'+str(neighbor)
					neighbor_rows.append(cr)
					neighbor_targets.append(neighbor)

		elif (self.target).startswith('Kepler') or (self.target).startswith('kepler') or (self.target).startswith("KEPLER"):
			for cr in check_rows:
				if (np.array(mast_data['kepler_name'])[cr][7:-2] == NEA_targetname[7:-2]) and (cr != rowidx):
					print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['kepler_name'][cr]))
					neighbor = str(mast_data['kepler_name'][cr])
					if ' ' in neighbor:
						neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
					neighbor_rows.append(cr)
					neighbor_targets.append(neighbor)

		elif (self.target).startswith('K2') or (self.target).startswith('k2'):
			for cr in check_rows:
				if (np.array(mast_data['pl_name'])[cr][3:-2] == NEA_targetname[3:-2]) and (np.array(mast_data['pl_name'])[cr] != NEA_targetname):
					print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['pl_name'][cr]))
					neighbor = str(mast_data['pl_name'][cr])
					if ' ' in neighbor:
						neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
					neighbor_rows.append(cr)
					neighbor_targets.append(neighbor)		

		else:
			if (self.telescope == "TESS") or (self.telescope == "Tess") or (self.telescope == 'tess'):
				if self.target.startswith("TIC"):
					for cr in check_rows:
						if (str(np.array(exofop_data['TIC ID'])[cr])[3:] == str(NEA_targetname)[3:]) and (cr != rowidx):
							print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(exofop_data['TIC ID'][cr]))
							neighbor = str(exofop_data['TIC ID'][cr])
							if ' ' in neighbor:
								neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
							neighbor_rows.append(cr)
							neighbor_targets.append('TIC '+str(neighbor))

				elif self.target.startswith('TOI'):
					for cr in check_rows:
						if (str(np.array(exofop_data['TOI'])[cr])[:-2] == str(NEA_targetname)[4:-2]) and (cr != rowidx):
							print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(exofop_data['TOI'][cr]))
							neighbor = str(exofop_data['TOI'][cr])
							if ' ' in neighbor:
								neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
							neighbor_rows.append(cr)
							neighbor_targets.append('TOI-'+str(neighbor))						

				else:
					for cr in check_rows:
						if (np.array(mast_data['pl_name'])[cr][3:-2] == NEA_targetname[3:-2]) and (np.array(mast_data['pl_name'])[cr] != NEA_targetname):
							print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(mast_data['pl_name'][cr]))
							neighbor = str(mast_data['pl_name'][cr])
							if ' ' in neighbor:
								neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
							neighbor_rows.append(cr)
							neighbor_targets.append(neighbor)									

		self.neighbors = neighbor_targets



	def get_neighbors(self, clobber_lc='y', save_to_file='y'):
		### this function will download grab all the information about your target's neighbors.
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
					if (self.telescope == 'kepler') or (self.telescope == 'Kepler') or (self.telescope == 'KEPLER'):
						if neighbor.startswith('Kepler') or neighbor.startswith('kepler'):
							neighbor_key = 'k'+str(neighbor[7:])
						else:
							### for now this is just a KOI!
							neighbor_key = neighbor
							while neighbor_key.startswith('K') or neighbor_key.startswith('0'):
								neighbor_key = neighbor_key[1:]
							neighbor_key = 'k'+str(neighbor_key)
					elif (self.telescope == 'k2') or (self.telescope == 'K2'):
						neighbor_key = 'k2_'+str(neighbor_key[3:])

					elif (self.telescope == 'tess') or (self.telescope == 'TESS') or (self.telescope == 'Tess'):
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
			lcfile = open(savepath+'/'+self.target+"_lightcurve.tsv", mode='r')
			lcfile_new = open(savepath+"/"+self.target+"_lc_temp.tsv", mode='w')

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
			os.system('mv '+savepath+'/'+self.target+'_lc_temp.tsv '+savepath+'/'+self.target+'_lightcurve.tsv')

		self.all_transit_times = np.array(neighbor_transit_times)
		self.all_transit_list = neighbor_transit_list
		self.all_transit_IDs = neighbor_transit_ID 



	def prep_for_CNN(self, save_lc='y', window=6, cnn_len=493, exclude_neighbors='y', flag_neighbors='y', show_plot='n'):
		### this function will produce an arrayy that's ready to be fed to a CNN for moon finding. 
		cnnlc_path = moonpydir+'/cnn_lcs'
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

			localpath = cnnlc_path+'/'+self.target+"_transit"+str(taunum)+'_cnnstack.npy'
			localpath_list.append(localpath)

			np.save(localpath, cnn_stack)

			if show_plot == 'y':
				plt.scatter(cnn_stack[0], cnn_stack[3], facecolor='LightCoral', edgecolor='k', s=10)
				if flag_neighbors == 'y':
					neighbor_transit_idxs = np.where(flag_array == 1)[0]
					plt.scatter(cnn_stack[0][neighbor_transit_idxs], cnn_stack[3][neighbor_transit_idxs], c='g', marker='x', s=15)
				plt.xlabel("BKJD")
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


	def plot_lc(self, facecolor='LightCoral', edgecolor='k', errorbar='n', quarters='all', folded='n', include_flagged='n', detrended='y', show_errors='n', show_neighbors='n'):
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
					plt.scatter(stitched_times[neighbor_transit_idxs], stitched_fluxes[neighbor_transit_idxs], s=10, marker='x')					


		elif folded == 'y':
			try:
				self.fold(detrended=detrended)
			except:
				self.get_properties(locate_neighbor='n')
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


