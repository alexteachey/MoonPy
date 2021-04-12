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
moonpydir = moonpydir[:moonpydir.find('/_mp_attributes.py')]



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
			



def find_planet_row(self, row_known='n'):
	download_new = 'n'
	if row_known == 'y':
		mast_rowidx = self.mast_rowidx
		if self.telescope.lower() != "kepler":
			exofop_rowidx = self.exofop_rowidx


	#### NEW CODE -- stolen from "planet_list_and_data_retriever.py"
	### check the file created times:
	kep_mast_address = moonpydir+'/cumkois_mast.txt'
	if os.path.exists(kep_mast_address):
		kep_mast_fct = os.path.getctime(kep_mast_address) ### file created time
	else:
		kep_mast_fct = 0
	kep_fop_address = moonpydir+'/kepler_exofop_targets.csv'
	if os.path.exists(kep_fop_address):
		kep_fop_fct = os.path.getctime(kep_fop_address) ### file created time
	else:
		kep_fop_fct = 0

	k2_mast_address = moonpydir+'/cumk2ois_mast.txt'
	if os.path.exists(kep_mast_address):
		k2_mast_fct = os.path.getctime(kep_mast_address) ### file created time
	else:
		k2_mast_fct = 0
	k2_fop_address = moonpydir+'/k2_exofop_targets.csv'
	if os.path.exists(k2_fop_address):
		k2_fop_fct = os.path.getctime(kep_fop_address) ### file created time
	else:
		k2_fop_fct = 0

	tess_mast_address = moonpydir+'/mast_cumtois.txt'
	if os.path.exists(tess_mast_address):
		tess_mast_fct = os.path.getctime(tess_mast_address) ### file created time
	else:
		tess_mast_fct = 0
	tess_fop_address = moonpydir+'/tess_exofop_targets.csv'
	if os.path.exists(tess_fop_address):
		tess_fop_fct = os.path.getctime(tess_fop_address) ### file created time
	else:
		tess_fop_fct = 0


	current_time = time.time()

	#### KEPLER HANDLING
	if current_time - kep_mast_fct > 86400: ### one day old
		print("DOWNLOADING Kepler MAST file (once per day)...")
		os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,koi_smass,koi_smass_err1,koi_smass_err2,koi_model_snr,ra,dec&order=kepoi_name&format=ascii" -O "'+kep_mast_address+'"')

		print(" ")
	if current_time - kep_fop_fct > 86400:
		print("DOWNLOADING Kepler ExoFOP file (once per day)...")
		os.system('wget --tries=1 --user=teachey --password=Tipiu2ExoFOP "https://exofop.ipac.caltech.edu/kepler/download_summary_csv.php?sort=koi" -O "'+kep_fop_address+'"')

		print(' ')

	#### K2 HANDLING
	if current_time - k2_mast_fct > 86400:
		print('DOWNLOADING K2 MAST file (once per day)...')
		os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=k2candidates&select=epic_name,epic_candname,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandep,pl_trandeperr1,pl_trandeperr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_eqt,pl_eqterr1,pl_eqterr2,pl_rade,pl_radeerr1,pl_radeerr2,st_rad,st_raderr1,st_raderr2,ra,dec&order=epic_name&format=ascii" -O "'+k2_mast_address+'"')
		print(' ')
	if current_time - k2_fop_fct > 86400:	
		print("DOWNLOADING K2 ExoFOP file (once per day)...")
		os.system('wget --tries=1 "https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All&sort=target" -O "'+k2_fop_address+'"')
		#os.system('wget --tries=1 '+k2_fop_address+' "https://exofop.ipac.caltech.edu/k2/download_summary_csv/php?camp=All&sort=target"')
		print(' ')


	#### TESS HANDLING
	if current_time - tess_mast_fct > 86400:
		print('DOWNLOADING TESS MAST file (once per day)...')
		os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=pl_hostname,pl_letter,pl_name,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_imppar,pl_impparerr1,pl_impparerr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_rade,pl_radeerr1,pl_radeerr2,ra,dec&order=pl_hostname&format=ascii" -O "'+tess_mast_address+'"')			
		print(' ')
	if current_time - tess_fop_fct > 86400:
		print('DOWNLOADING TESS ExoFOP file (once per day)...')
		os.system('wget --tries=1 "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv" -O "'+tess_fop_address+'"')				
		print(' ')


	### USER INPUT HANDLING
	if self.telescope.lower() == 'user':
		mast_data = ascii.read(kep_mast_address)
		mast_columns = mast_data.columns
		exofop_data = pandas.read_csv(kep_fop_address, header=18)
		exofop_columns = exofop_data.columns
		row_known = 'y'
		self.mast_rowidx = np.nan 
		self.exofop_rowidx = np.nan 

	### KEPLER HANDLING 
	elif (self.telescope.lower() == 'kepler'):
		mast_data = ascii.read(kep_mast_address)
		mast_columns = mast_data.columns
		try:
			exofop_data = pandas.read_csv(kep_fop_address, header=18)
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
				mast_data = ascii.read(kep_mast_address) ### overwrite before! won't be found in cumkois!
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
		mast_data = ascii.read(k2_mast_address)
		mast_columns = mast_data.columns
		try:
			exofop_data = pandas.read_csv(k2_fop_address, header=10)
			exofop_columns = exofop_data.columns
		except:
			exofop_data = None
			exofop_columns = np.array([])
			print('K2 ExoFOP file not loadable. Possibly corrupted. ')

		if row_known == 'n':
			if str(self.target).startswith('K2-'):
				target_letter = str(self.target[-1])
				if ' ' in self.target:
					NEA_targetname = self.target
				else:
					NEA_targetname = str(self.target[:-1])+' '+str(self.target[-1])
				mast_rowidx = np.where(mast_data['pl_name'] == NEA_targetname)[0]
				exofop_rowidx = np.nan
				print('mast_rowidx = ', mast_rowidx)
				print('exofop_rowidx = ', exofop_rowidx) 

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
		mast_data = ascii.read(tess_mast_address)
		mast_columns = mast_data.columns
		exofop_data = pandas.read_csv(tess_fop_address)
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
		target_smass, target_smass_uperr, target_smass_lowerr = mast_data['koi_smass'][mast_rowidx], mast_data['koi_smass_err1'][mast_rowidx], mast_data['koi_smass_err2'][mast_rowidx]
		#target_snr = mast_data['koi_model_snr']


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
			target_tau0, target_tau0_uperr, target_tau0_lowerr = np.array(exofop_data['Epoch (BJD)'])[exofop_rowidx], np.array(exofop_data['Epoch (BJD) err'])[exofop_rowidx], np.array(exofop_data['Epoch (BJD) err'])[exofop_rowidx]

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


	if (self.telescope.lower() == 'k2') and (float(target_tau0) > 2454833):
		try:
			self.tau0 = float(target_tau0) - 2454833 #### convert to BKJD
		except:
			self.tau0 = np.nan


	elif (self.telescope.lower() == 'tess') and (float(target_tau0) > 2454833):
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
		#self.snr = float(target_snr)

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

	try:
		self.smass = float(target_smass)
		self.smass_err = (float(target_smass_lowerr), float(target_smass_uperr))
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
		#continue_query = input('Do you wish to continue? y/n: ')
		#if continue_query != 'y':
		#	raise Exception('you opted not to continue.')
		while (transit_midtimes[-1] - self.period) > np.nanmin(np.hstack(self.times)):
			print('appending transit_midtime: ', transit_midtimes[-1] - self.period)
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
		print('mast_rowidx, NEA_targetname = ', mast_rowidx, NEA_targetname)

		if len(mast_rowidx) == 0:
			### means this object is not in the cumulative KOI list.
			print('It appears this object is not in the cumulative KOI list.')
			print('checking all '+str(len(mast_data['kepoi_name']))+' MAST rows...')
			#### you're going to need to check all mast rows!
			check_mast_rows = np.arange(0,len(mast_data['kepoi_name']),1)
		else:
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
				target_term = self.NEA_targetname[:-1]
				search_term = np.array(mast_data['kepoi_name'])[cr][:-1]
				if (search_term == target_term) and (cr != mast_rowidx):
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
				target_term = self.target[7:-1]
				search_term = np.array(mast_data['kepler_name'])[cr][7:-2]
				if ((len(mast_rowidx) == 0) and (search_term == target_term)) or ((search_term == target_term) and (cr != mast_rowidx)):
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
				target_term = NEA_targetname[3:-2]
				search_term = np.array(mast_data['pl_name'])[cr][3:-2]
				if (search_term == target_term) and (np.array(mast_data['pl_name'])[cr] != NEA_targetname):
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
						neighbor_key = 'K'+str(neighbor[7:])
					else:
						### for now this is just a KOI!
						neighbor_key = neighbor
						while neighbor_key.startswith('K') or neighbor_key.startswith('0'):
							neighbor_key = neighbor_key[1:]
						neighbor_key = 'K'+str(neighbor_key)
				elif (self.telescope.lower() == 'k2'):
					neighbor_key = 'K2_'+str(neighbor_key[3:])

				elif (self.telescope.lower() == 'tess'):
					neighbor_key = neighbor 

				### now generate the LC object -- note that this will download duplicate light curves!
				### for now, just delete them after downloading... should come up with a better solution.
				print('calling MoonpyLC for neighbor: ', neighbor)

				neighbor_dict[neighbor_key] = MoonpyLC(targetID=neighbor, is_neighbor=is_neighbor, download='y', clobber='n')

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
			lcfile = open(self.savepath+'/'+self.target+"_"+self.telescope+"_lightcurve.tsv", mode='r')
		except:
			lcfile = open(self.savepath+'/'+self.target+'_lightcurve.tsv', mode='r') ### for older files.
		lcfile_new = open(self.savepath+"/"+self.target+"_"+self.telescope+"_lc_temp.tsv", mode='w')

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
		#os.system('mv '+self.savepath+'/'+self.target+'_lc_temp.tsv '+savepath+'/'+self.target+'_lightcurve.tsv')
		os.system('mv '+self.savepath+'/'+self.target+'_'+self.telescope+'_lc_temp.tsv '+self.savepath+'/'+self.target+'_'+self.telescope+'_lightcurve.tsv')

	self.neighbor_transit_times = np.array(neighbor_transit_times)
	self.neighbor_transit_list = neighbor_transit_list
	self.neighbor_transit_IDs = neighbor_transit_ID 






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




def get_future_transits(self, num_transits=20, output_format='datetime', native_format=None):
	#### uses last tau to project future transit times
	#### format can be 'BJD', 'BKJD', 'BTJD', or 'datetime'
	last_tau = self.taus[-1]
	print('producing '+str(num_transits)+' future transits. Barycentric corrections not considered.')

	if self.telescope.lower() == 'kepler':
		native_format='BKJD'
	elif self.telescope.lower() == 'tess':
		native_format=='BTJD'

	if native_format.lower() == 'bkjd':
		last_tau_BJD = last_tau + 2454833.0
	elif native_format.lower() == 'btjd':
		last_tau_BJD = last_tau + 2457000.0
	elif (native_format.lower() == 'bjd') or (native_format.lower() == 'jd'):
		last_tau_BJD = last_tau 

	future_taus = []
	for i in np.arange(1,num_transits+1,1):
		future_taus.append(last_tau_BJD + (i*self.period))
	future_taus = np.array(future_taus)
	future_taus = Time(future_taus, format='jd', scale='utc')

	if output_format.lower() == 'datetime':
		output_taus = future_taus.isot 
	elif (output_format.lower() == 'bjd') or (output_format.lower() == 'jd'):
		output_taus = future_taus.jd 
	elif (output_format.lower() == 'bkjd') or (output_format.lower() == 'kjd'):
		output_taus = future_taus.jd - 2454833
	elif (output_format.lower() == 'btjd') or (output_format.lower() == 'tjd'):
		output_taus = future_taus.jd - 2457000

	
	### return future_taus in both native format, and the requested format
	#for opt, ft in zip(output_taus, future_taus):
	#		print(str(opt)+' --- '+str(ft))

	return future_taus, output_taus















