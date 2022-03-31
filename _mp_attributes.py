from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
from astropy.io import ascii
from astropy.time import Time
import time
import pandas
import traceback
from astroquery.simbad import Simbad 
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
from astropy.timeseries import LombScargle
import socket 

#### BELOW ARE MOONPY PACKAGES
from moonpy import *
from mp_tools import *
from mp_lcfind import *
from mp_detrend import untrendy_detrend, cofiam_detrend, george_detrend, medfilt_detrend, polyAM_detrend
from mp_batman import run_batman
from mp_fit import mp_multinest, mp_emcee
from mp_plotter import *
from cofiam import max_order
from pyluna import prepare_files
from mp_tpf_examiner import *
from scipy.interpolate import interp1d 
from moonpy import *

plt.rcParams["font.family"] = 'serif'

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/_mp_attributes.py')]


#### test the moonpy import right here
#print("TESTING THE MOONPY IMPORT (_mp_attributes.py)")
#dummy = MoonpyLC(targetID='Kepler-1625b')


def find_transit_quarters(self, locate_neighbor='n'):
	print('calling _mp_attributes.py/find_transit_quarters().')
	self.get_properties(locate_neighbor=locate_neighbor) ### calls find_planet_row() within.
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

			quarter_times = np.array(quarter_times)

			quarter_transit = 'n'
			for tau in self.taus:
				if (tau >= np.nanmin(quarter_times)) and (tau <= np.nanmax(quarter_times)):
					### there's a transit in this quarter!
					quarter_transit = 'y'
					break

			quarter_transit_dict[quarter] = quarter_transit

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
	print(self.quarter_transit_dict)
	print(' ')


def find_aliases(self):
	print('calling _mp_attributes.py/find_aliases().')
	target_aliases = []
	try:
		alias_search_results = Simbad.query_objectids(self.target)
		for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
			target_alias = alias_search_results[alidx][0]
			target_aliases.append(target_alias)
		if len(target_aliases) == 0:
			print('NO SIMBAD ALIASES FOUND. RETURNING TARGET ONLY.')
			target_aliases = [self.target]
		self.aliases = np.array(target_aliases)
		print(self.aliases)
		print( )
	except:
		print('UNABLE TO IDENTIFY ALIASES.')




def get_coords(self):
	print('calling _mp_attributes.py/get_coords().')
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
					self.RA = self.NEA_data['ra'][self.NEA_rowidx]
					self.Dec = self.NEA_data['dec'][self.NEA_rowidx]
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
					self.RA = np.array(self.exofop_data['RA'][self.exofop_rowidx])[0]
					self.Dec = np.array(self.exofop_data["DEC"][self.exofop_rowidx])[0]
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
							self.RA = np.array(self.exofop_data['RA'][self.exofop_rowidx])[0]
							self.Dec = np.array(self.exofop_data['Dec'][self.exofop_rowidx])[0]
							self.ticnum = np.array(self.exofop_data['TIC ID'][self.exofop_rowidx])[0]
						except:
							pass

				elif self.target_type.lower() == 'tic':
					try:
						simbad_query = Simbad.query_object('TIC '+self.target)[0] 
						target_name = 'TIC '+str(self.target) ### only update if this worked!
						self.RA = simbad_query['RA']
						self.Dec = simbad_query['DEC']
						self.ticnum = target_name
						#target_in_simbad = 'y'
					except:
						#target_in_simbad = 'n'
						try:
							self.RA = np.array(self.exofop_data['RA'][self.exofop_rowidx])[0]
							self.Dec = np.array(self.exofop_data['Dec'][self.exofop_rowidx])[0]
							self.ticnum = np.array(self.exofop_data['TIC ID'][self.exofop_rowidx])[0]
						except:
							pass





def get_exofop_credentials():

	### check for an exofop credentials file
	if os.path.exists(moonpydir+'/exofop_credentials.txt'):
		#### open and read it
		exofop_credentials = open(moonpydir+'/exofop_credentials.txt', mode='r')
		exofop_username, exofop_pw = exofop_credentials.readline().split()
		exofop_credentials.close()

	else:
		#### make the user input these, and save the file
		print(' ')
		print(' ')
		print('NOTE: MoonPy makes use of the ExoFOP databases. ')
		print('Users are advised to set up an account if they do not have one already. ')
		print(' ')
		print(' ')


		exofop_username = input('Please enter your exofop username (first time user only), or press ENTER if you do not have one: ')
		if len(exofop_username) > 0:
			exofop_pw = input('Please enter your exofop password: ')

			exofop_credentials = open(moonpydir+'/exofop_credentials.txt', mode='w')
			exofop_credentials.write(exofop_username+' '+exofop_pw)
			exofop_credentials.close()

	return exofop_username, exofop_pw 





#def get_databases(self):
def get_databases(target_prefix):

	current_time = time.time()


	#### KEPLER HANDLING
	#if (self.telescope.lower() == 'kepler') :

	#### NEW FUNCTIONALITY -- ONE NEA TABLE FOR KEPLER, K2 AND TESS. 
	kep_koi_address = moonpydir+'/cumkois_mast.csv'
	kep_NEA_address = moonpydir+'/confirmed_planets_mast.csv'
	
	if os.path.exists(kep_NEA_address):
		kep_NEA_fct = os.path.getctime(kep_NEA_address) ### file created time
	else:
		kep_NEA_fct = 0

	if os.path.exists(kep_koi_address):
		kep_koi_fct = os.path.getctime(kep_koi_address) ### file created time
	else:
		kep_koi_fct = 0


	if current_time - kep_NEA_fct > 86400: ### one day old
		#### CHECK WITH USER ABOUT RE-DOWNLOADING.
		download_new_csvs = input('Planet databases are missing or more than one day out of date. Do you want to download the new versions? (Takes a few minutes). y/n: ')
	else:
		download_new_csvs = 'n'


	if download_new_csvs == 'y':

		exofop_username, exofop_pw = get_exofop_credentials()


		print("DOWNLOADING Kepler MAST file (once per day)...")

		desired_columns = np.array(['pl_name', 'pl_letter', 'hostname', 'hd_name', 'hip_name',
	       'tic_id', 'disc_pubdate', 'disc_year', 'discoverymethod', 'disc_locale', 'disc_facility', 'disc_telescope',
	       'ra', 'dec', 'glon', 'glat', 'elon', 'elat', 
	       'pl_orbper', 'pl_orbpererr1', 'pl_orbpererr2', 'pl_orbperlim', ### orbital period
	       'pl_orblpererr1','pl_orblper', 'pl_orblpererr2', 'pl_orblperlim', #### argument of periastron
	       'pl_orbsmax', 'pl_orbsmaxerr1', 'pl_orbsmaxerr2', 'pl_orbsmaxlim', #### semimajor axis
	       'pl_orbincl', 'pl_orbinclerr1', 'pl_orbinclerr2', 'pl_orbincllim', #### orbital inclination
	       'pl_orbtper', 'pl_orbtpererr1', 'pl_orbtpererr2', 'pl_orbtperlim', #### epoch of periastron
	       'pl_orbeccen', 'pl_orbeccenerr1', 'pl_orbeccenerr2', 'pl_orbeccenlim', #### eccentricity
	       'pl_eqt', 'pl_eqterr1', 'pl_eqterr2', 'pl_eqtlim', #### equilibrium temperature
	       'pl_occdep', 'pl_occdeperr1', 'pl_occdeperr2', 'pl_occdeplim',  ### occultation depth
	       'pl_insol', 'pl_insolerr1', 'pl_insolerr2', 'pl_insollim', ### insolation
	       'pl_dens', 'pl_denserr1', 'pl_denserr2', 'pl_denslim', #### density 
	       'pl_trandep', 'pl_trandeperr1', 'pl_trandeperr2', 'pl_trandeplim', ### transit depth
	       'pl_tranmid', 'pl_tranmiderr1', 'pl_tranmiderr2', 'pl_tranmidlim', 'pl_tsystemref', #### transit midtime
	       'sy_pmdec', 'sy_pmdecerr1', 'sy_pmdecerr2', 
	       'sy_plx', 'sy_plxerr1', 'sy_plxerr2', 'sy_dist', 'sy_disterr1', 'sy_disterr2', 
	       'sy_bmag', 'sy_bmagerr1', 'sy_bmagerr2', 'sy_vmag', 'sy_vmagerr1','sy_vmagerr2', 
	       'sy_jmag', 'sy_jmagerr1', 'sy_jmagerr2', 'sy_hmag','sy_hmagerr1', 'sy_hmagerr2', 
	       'sy_kmag', 'sy_kmagerr1','sy_kmagerr2', 'sy_umag', 'sy_umagerr1', 'sy_umagerr2', 
	       'sy_rmag', 'sy_rmagerr1', 'sy_rmagerr2', 'sy_imag', 'sy_imagerr1',
	       'sy_imagerr2', 'sy_zmag', 'sy_zmagerr1', 'sy_zmagerr2', 
	       'sy_gaiamag', 'sy_gaiamagerr1','sy_gaiamagerr2', 'sy_tmag', 'sy_tmagerr1', 'sy_tmagerr2',
	       'pl_controv_flag',
	       'st_metratio', 'st_spectype', 'sy_kepmag', 'sy_kepmagerr1',
	       'sy_kepmagerr2', 'st_rotp', 'st_rotperr1', 'st_rotperr2',
	       'st_rotplim', 'pl_projobliq', 'pl_projobliqerr1',
	       'pl_projobliqerr2', 'pl_projobliqlim', 'gaia_id', 'cb_flag',
	       'pl_trandur', 'pl_trandurerr1', 'pl_trandurerr2', 'pl_trandurlim',
	       'pl_rvamp', 'pl_rvamperr1', 'pl_rvamperr2', 'pl_rvamplim',
	       'pl_radj', 'pl_radjerr1', 'pl_radjerr2', 'pl_radjlim', 'pl_rade',
	       'pl_radeerr1', 'pl_radeerr2', 'pl_radelim', 'pl_ratror',
	       'pl_ratrorerr1', 'pl_ratrorerr2', 'pl_ratrorlim', 'pl_ratdor',
	       'pl_trueobliq', 'pl_trueobliqerr1', 'pl_trueobliqerr2',
	       'pl_trueobliqlim', 'sy_icmag', 'sy_icmagerr1', 'sy_icmagerr2',
	       'dkin_flag', 'pl_ratdorerr1', 'pl_ratdorerr2', 'pl_ratdorlim',
	       'pl_imppar', 'pl_impparerr1', 'pl_impparerr2', 'pl_impparlim',
	       'pl_bmassj', 'pl_bmassjerr1', 'pl_bmassjerr2', 'pl_bmassjlim',
	       'pl_bmasse', 'pl_bmasseerr1', 'pl_bmasseerr2', 'pl_bmasselim',
	       'pl_bmassprov', 'st_teff', 'st_tefferr1', 'st_tefferr2',
	       'st_tefflim', 'st_met', 'st_meterr1', 'st_meterr2', 'st_metlim',
	       'st_radv', 'st_radverr1', 'st_radverr2', 'st_radvlim', 'st_vsin',
	       'st_vsinerr1', 'st_vsinerr2', 'st_vsinlim', 'st_lum', 'st_lumerr1',
	       'st_lumerr2', 'st_lumlim', 'st_logg', 'st_loggerr1', 'st_loggerr2',
	       'st_logglim', 'st_age', 'st_ageerr1', 'st_ageerr2', 'st_agelim',
	       'st_mass', 'st_masserr1', 'st_masserr2', 'st_masslim', 'st_dens',
	       'st_denserr1', 'st_denserr2', 'st_denslim', 'st_rad', 'st_raderr1',
	       'st_raderr2', 'st_radlim', 'ttv_flag', 'ptv_flag', 'tran_flag',
	       'rv_flag', 'ast_flag', 'obm_flag', 'micro_flag', 'etv_flag',
	       'ima_flag', 'pul_flag', 'sy_snum', 'sy_pnum', 'sy_mnum',
	       'st_nphot', 'st_nrvc', 'st_nspec', 'pl_nespec', 'pl_ntranspec',
	       'pl_nnotes', 'sy_pm', 'sy_pmerr1', 'sy_pmerr2', 'sy_pmra',
	       'sy_pmraerr1', 'sy_pmraerr2', 'x', 'y', 'z', 'htm20'],
	      dtype='<U20')

		print('number of desired columns: ', len(desired_columns))
		
		desired_cols = ''
		for ndc,dc in enumerate(desired_columns):

			desired_cols = desired_cols+dc+','
		if desired_cols[-1] == ',':
			desired_cols = desired_cols[:-1]

		print(' ')
		print(' ')
		print('desired_cols string: ')
		print(desired_cols)
		print(' ')
		print(' ')

		#NEW_query =  '"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+'+desired_cols+'+from+pscomppars&format=csv" -O "'+kep_NEA_address+'"'
		NEW_query =  '"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+pscomppars+order+by+pl_name+asc&format=csv" -O "'+kep_NEA_address+'"'

		#os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,koi_smass,koi_smass_err1,koi_smass_err2,koi_model_snr,ra,dec&order=kepoi_name&format=ascii" -O "'+kep_NEA_address+'"')

		os.system('wget -v --tries=1 '+NEW_query)
		print(" ")




	if (current_time - kep_koi_fct > 86400) and (download_new_csvs == 'y'): ### one day old
		print("DOWNLOADING KOI MAST file (once per day)...")
		OLD_query = 'wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,koi_smass,koi_smass_err1,koi_smass_err2,koi_model_snr,ra,dec&order=kepoi_name&format=ascii" -O "'+kep_koi_address+'"'
		os.system('wget -v --tries=1 '+OLD_query)	




	#### LOAD THE DATA IN
	confirmed_NEA_data = ascii.read(kep_NEA_address)
	#NEA_data = pandas.read_csv(kep_NEA_address)
	confirmed_NEA_columns = confirmed_NEA_data.columns
	koi_NEA_data = ascii.read(kep_koi_address)
	#NEA_data = pandas.read_csv(kep_NEA_address)
	koi_NEA_columns = koi_NEA_data.columns

	#if self.target.lower().startswith('kepler'):
	if target_prefix.lower() == 'kepler':
		NEA_data = confirmed_NEA_data
		NEA_columns = confirmed_NEA_columns

	elif (target_prefix.lower() == 'koi') or (target_prefix.lower() == 'kic'):
		NEA_data = koi_NEA_data
		NEA_columns = koi_NEA_columns 

	elif (target_prefix.lower() == 'k2') or (target_prefix.lower() == 'epic'):
		NEA_data = confirmed_NEA_data
		NEA_columns = confirmed_NEA_columns 

	elif (target_prefix.lower() == 'toi') or (target_prefix.lower() == 'tic') or (target_prefix.lower() == 'tess'):
		NEA_data = confirmed_NEA_data 
		NEA_columns = confirmed_NEA_columns




	if (target_prefix.lower() == 'kepler') or (target_prefix.lower() == 'koi') or (target_prefix.lower() == 'kic'):
		kep_fop_address = moonpydir+'/kepler_exofop_targets.csv'
	
		if os.path.exists(kep_fop_address):
			kep_fop_fct = os.path.getctime(kep_fop_address) ### file created time
		else:
			kep_fop_fct = 0

		if (current_time - kep_fop_fct > 86400) and (download_new_csvs == 'y'):
			exofop_username, exofop_pw = get_exofop_credentials()

			print("DOWNLOADING Kepler ExoFOP file (once per day)...")
			os.system('wget --tries=1 --user='+exofop_username+' --password='+exofop_pw+' "https://exofop.ipac.caltech.edu/kepler/download_summary_csv.php?sort=koi" -O "'+kep_fop_address+'"')	
			print(' ')

		try:
			exofop_data = pandas.read_csv(kep_fop_address, header=18)
			exofop_columns = exofop_data.columns
		except:
			exofop_data = None
			exofop_columns = np.array([])			
			print('EXCEPTION ENCOUNTERED... could not load exofop_data.')


	elif (target_prefix.lower() == 'k2') or (target_prefix.lower() == 'epic'):

		k2_fop_address = moonpydir+'/k2_exofop_targets.csv'
		if os.path.exists(k2_fop_address):
			k2_fop_fct = os.path.getctime(k2_fop_address) ### file created time
		else:
			k2_fop_fct = 0
			#download_new_csvs = 'y' 


		if (current_time - k2_fop_fct > 86400):	
			download_new_csvs = input('Planet databases are missing or more than one day out of date. Do you want to download the new versions? (Takes a few minutes). y/n: ')


			if download_new_csvs == 'y':

				exofop_username, exofop_pw = get_exofop_credentials()	

				print("DOWNLOADING K2 ExoFOP file (once per day)...")
				os.system('wget --tries=1 --user='+exofop_username+' --password='+exofop_pw+' "https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All&sort=target" -O "'+k2_fop_address+'"')
				#os.system('wget --tries=1 --user='+exofop_username+' --password='+exofop_pw+' "https://exofop.ipac.caltech.edu/kepler/download_summary_csv.php?sort=koi" -O "'+kep_fop_address+'"')	

				print(' ')

		try:
			exofop_data = pandas.read_csv(k2_fop_address, header=10)
			exofop_columns = exofop_data.columns
		except:
			exofop_data = None
			exofop_columns = np.array([])
			print('K2 ExoFOP file not loadable. Possibly corrupted. ')



	elif (target_prefix.lower() == 'tic') or (target_prefix.lower() == 'toi') or (target_prefix.lower() == 'tess'):

		tess_fop_address = moonpydir+'/tess_exofop_targets.csv'
		if os.path.exists(tess_fop_address):
			tess_fop_fct = os.path.getctime(tess_fop_address) ### file created time
		else:
			tess_fop_fct = 0

		if (current_time - tess_fop_fct > 86400):

			download_new_csvs = input('Planet databases are missing or more than one day out of date. Do you want to download the new versions? (Takes a few minutes). y/n: ')

			if download_new_csvs == 'y':

				exofop_username, exofop_pw = get_exofop_credentials()

				print('DOWNLOADING TESS ExoFOP file (once per day)...')
				#os.system('wget --tries=1 "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv" -O "'+tess_fop_address+'"')			
				os.system('wget --tries=1 --user='+exofop_username+' --password='+exofop_pw+' "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv" -O "'+tess_fop_address+'"')						
				print(' ')

		#NEA_data = ascii.read(NEA_address)
		#NEA_columns = NEA_data.columns 
		exofop_data = pandas.read_csv(tess_fop_address)
		exofop_columns = exofop_data.columns


	return NEA_data, NEA_columns, exofop_data, exofop_columns 











def find_planet_row(self, alias=None, row_known='n'):

	print('calling _mp_attributes.py/find_planet_row().')
	if alias != None:
		row_target = alias
	else:
		row_target = self.target 

	download_new = 'n'

	#### stolen from "planet_list_and_data_retriever.py"
	### check the file created times:

	if self.telescope.lower() == 'kepler':
		kep_koi_address = moonpydir+'/cumkois_mast.csv'
		kep_NEA_address = moonpydir+'/confirmed_planets_mast.csv'

		if row_known == 'n':
			if str(row_target).lower().startswith('kepler'):

				target_letter = str(row_target)[-1]
				if ' ' in row_target: ### already in the correct format, with a space between the letter.
					NEA_targetname = row_target
				else: #### there isn't a space, so add it!
					NEA_targetname = row_target[:-1]+' '+target_letter
				#NEA_rowidx = np.where(np.char.lower(NEA_data['kepler_name']) == NEA_targetname.lower())[0]
				if "NEA_rowidx" in dir(self):
					NEA_rowidx = self.NEA_rowidx
				else:
					NEA_rowidx = np.where(np.char.lower(self.NEA_data['pl_name']) == NEA_targetname.lower())[0]		
				print("NEA_rowidx = ", NEA_rowidx)		


			elif str(row_target).lower().startswith('kic'):
				NEA_targetname = int(row_target[4:])
				if 'NEA_rowidx' in dir(self):
					NEA_rowidx = self.NEA_rowidx 
				else:
					#### THIS HAS TO CHANGE!
					NEA_rowidx = np.where(self.NEA_data['kepid'] == NEA_targetname)[0]

			elif str(row_target).lower().startswith('koi'):
				NEA_targetname = str(row_target[4:])
				if len(NEA_targetname) == 7: ### of the form 5084.01
					NEA_targetname = 'K0'+str(NEA_targetname)
				elif len(NEA_targetname) == 6: ### of the form 163.01
					NEA_targetname = 'K00'+str(NEA_targetname)
				elif len(NEA_targetname) == 5: ### of the form 23.01
					NEA_targetname = 'K000'+str(NEA_targetname)
				elif len(NEA_targetname) == 4: ### of the form 1.01
					NEA_targetname = 'K0000'+str(NEA_targetname)

				if "NEA_rowidx" in dir(self):
					NEA_rowidx = self.NEA_rowidx
				else:
					NEA_rowidx = np.where(np.char.lower(self.NEA_data['kepoi_name']) == NEA_targetname.lower())[0]
				self.NEA_targetname = NEA_targetname


			else: ### was observed by Kepler but you don't know it's KOI/KIC or Kepler name:
				try:
					float(row_target[-1]) ### if this works, query the NEA_data['pl_hostname'] because you end with a number!
					NEA_targetname = row_target 
					if "NEA_rowidx" in dir(self):
						NEA_rowidx = self.NEA_rowidx
					else:
						NEA_rowidx = np.where(np.char.lower(self.NEA_data['pl_hostname']) == row_target.lower())[0]

				except:
					### implies the last thing is a number
					target_letter = row_target[-1]
					if ' ' in row_target:
						NEA_targetname = row_target
					else:
						NEA_targetname = row_target[:-1]+' '+target_letter
					if "NEA_rowidx" in dir(self):
						NEA_rowidx = self.NEA_rowidx			
					else:			
						NEA_rowidx = np.where(np.char.lower(self.NEA_data['pl_name']) == NEA_targetname.lower())[0]


		elif row_known == 'y':
			try:
				NEA_rowidx = self.NEA_rowidx
			except:
				pass
			try:
				exofop_rowidx = self.exofop_rowidx
			except:
				pass

		try:
			#self.exofop_rowidx = int(exofop_rowidx)
			exofop_rowidx = int(exofop_rowidx)
		except:
			exofop_rowidx = np.nan 
			#self.exofop_rowidx = exofop_rowidx

		try:
			NEA_rowidx = int(NEA_rowidx)
		except:
			NEA_rowidx = np.nan 
			#self.NEA_rowidx = np.nan 


		try:
			#return NEA_rowidx, NEA_data, NEA_targetname
			return NEA_rowidx, NEA_targetname
		except:
			#return NEA_rowidx, NEA_data, row_target
			return NEA_rowidx, row_target 


	
	#### K2 HANDLING
	elif self.telescope.lower() == 'k2':
		k2_NEA_address = moonpydir+'/cumk2ois_mast.txt'

		if row_known == 'n':
			if str(row_target).lower().startswith('k2'):
				target_letter = str(row_target[-1])
				if ' ' in row_target:
					NEA_targetname = row_target
				else:
					NEA_targetname = str(row_target[:-1])+' '+str(target_letter)

				if "NEA_rowidx" in dir(self):
					NEA_rowidx = self.NEA_rowidx
				else:	
					NEA_rowidx = np.where(np.char.lower(self.NEA_data['pl_name']) == NEA_targetname.lower())[0]

				exofop_rowidx = np.nan
				print('NEA_rowidx = ', NEA_rowidx)
				print('exofop_rowidx = ', exofop_rowidx) 
				print("number of rows matching this description = ", len(NEA_rowidx))

			elif str(row_target).lower().startswith('epic'):
				NEA_targetname = row_target[4:] ### just the numbers!
				if (NEA_targetname.startswith(' ')) or (NEA_targetname.startswith('-')):
					NEA_targetname = NEA_targetname[1:]
				
				try:
					if 'exofop_rowidx' in dir(self):
						exofop_rowidx = self.exofop_rowidx 
					else:
						exofop_rowidx = np.where(np.char.lower(self.exofop_data['EPIC ID']) == NEA_targetname.lower())[0]
				
				except:
					print('unable to extract the exofop_rowidx')
					exofop_rowidx = np.nan 

				try:
					if "NEA_rowidx" in dir(self):
						NEA_rowidx = self.NEA_rowidx	
					else:			
						NEA_rowidx = np.where(np.char.lower(self.NEA_data['epic_name']) == NEA_targetname.lower())[0]
				
				except:
					print('unable to extract the NEA_rowidx')
					NEA_rowidx = np.nan

			try:
				exofop_rowidx = int(exofop_rowidx)
				#self.exofop_rowidx = int(exofop_rowidx)
			except:
				exofop_rowidx = np.nan 
				#self.exofop_rowidx = exofop_rowidx

			try:
				#self.NEA_rowidx = int(NEA_rowidx)
				NEA_rowidx = int(NEA_rowidx)
			except:
				NEA_rowidx = np.nan
				#self.NEA_rowidx = NEA_rowidx 

		try:
			return NEA_rowidx, exofop_rowidx, NEA_targetname
		except:
			return NEA_rowidx, exofop_rowidx, row_target 





	#### TESS HANDLING
	elif self.telescope.lower() == 'tess':

		##### RESTRUCTURE -- if it starts with TOI, check if it has a number or letter at the e
		#toi_NEA_address = moonpydir+'/NEA_cumtois.txt' #### ALL THE TOIs that are NOT confirmed as planets. Of the form TOI-XXXX.01
		NEA_address = moonpydir+'/NEA_confirmed_planets.txt' #### confirmed planets (not limited to TESS)

		#if row_known == 'n':
			
		if str(row_target)[-1] in ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']: 
			#### WILL WANT TO LOOK IN THE NEA_data.
			target_letter = str(row_target[-1])
			### implies the last value is a letter.

			if ' ' in row_target:
				NEA_targetname = str(row_target) ### we want a space in here.
			else: ### put one in!
				NEA_targetname = str(row_target[:-1])+' '+str(target_letter)
		else: 
			#### will want to look in NEA_toi_data.
			NEA_targetname = str(row_target) #### will be something like TIC or TOI, maybe with a .01 on the end.


		if NEA_targetname.lower().startswith('tic'):
			### search the "TIC ID" column in ExoFOP 
			ticnum = NEA_targetname[3:] #### start the number after the 'TIC' prefix
			if (ticnum.startswith(' ')) or (ticnum.startswith('-')): #### make sure you've just got the number.
				ticnum = ticnum[1:] 

			if NEA_targetname[3] != ' ': ### make sure there's a space!
				NEA_targetname = NEA_targetname[:3]+' '+NEA_targetname[3:]

			#### QUERY EXOFOP
			print('looking for '+str(ticnum)+' in the exofop database.')	
			if 'exofop_rowidx' in dir(self):
				exofop_rowidx = self.exofop_rowidx
			else:				
				exofop_rowidx = np.where(np.array(self.exofop_data['TIC ID']).astype(int) == int(ticnum))[0]
			print('found '+str(len(exofop_idx)+' entries.'))
			
			if len(exofop_rowidx) > 1:
				exofop_rowidx = exofop_rowidx[0]
			print("exofop_rowidx = ", exofop_rowidx)

			#### IT'S NOT IN THE NEA database, don't even try.
			NEA_rowidx = np.nan		



		elif NEA_targetname.lower().startswith("toi"):
			toinum = NEA_targetname[3:] #### get rid of the intial 'TOI'
			if (toinum.startswith(' ')) or (toinum.startswith('-')): ### get rid of the blank or hyphen
				toinum = toinum[1:]

			print('looking for '+str(toinum)+' in the exofop database.')

			if str(toinum[-1]) in ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']:
				#### WILL BE FOUND IN THE NEA_data, but *not* in the exoFOP data.
				### search the TOI column 
				if NEA_targetname[-2] != ' ':
					#### add the space in there
					NEA_targetname[:-2]+' '+NEA_targetname[-2:]						

				try:
					if 'NEA_rowidx' in dir(self):
						NEA_rowidx = self.NEA_rowidx 
					else:
						NEA_rowidx = np.where(np.char.lower(np.array(self.NEA_data['pl_name'])) == NEA_targetname.lower())[0]
				except:
					print('TOI value not found in the NEA Confirmed Planets database.')
					NEA_rowidx = np.nan		
				print('NEA_rowidx = ', NEA_rowidx)	


			else:
				#### it's of the form TOI-XXXX.01, not TOI-XXXX b.
				try:
					if 'exofop_rowidx' in dir(self):
						exofop_rowidx = self.exofop_rowidx
					else:
						exofop_rowidx = np.where(np.array(self.exofop_data['TOI']).astype(str) == str(toinum))[0]
				except:
					try:
						if 'exofop_rowidx' in self.exofop_rowidx:
							exofop_rowidx = self.exofop_rowidx
						else:
							exofop_rowidx = np.where(np.array(self.exofop_data['TOI']).astype(str) == str(nospaces(toinum)))[0]
					except:
						print('TOI value not found in ExoFOP.')
						exofop_rowidx = np.nan
				print('exofop_rowidx = ', exofop_rowidx)



		else:
			### NOT A TIC or a TOI.... in this case, you have to look in the mast list! It's not listed in ExoFOP.
			#### have to find it in NEA_data
			if "NEA_rowidx" in dir(self):
				NEA_rowidx = self.NEA_rowidx 
			else:
				try:
					NEA_rowidx = np.where(np.char.lower(self.NEA_data['pl_name']) == NEA_targetname.lower())[0] ### in "confirmed" planets.
				except:
					NEA_rowidx = np.nan 

		try:
			#self.exofop_rowidx = int(exofop_rowidx)
			exofop_rowidx = int(exofop_rowidx)
		except:
			exofop_rowidx = np.nan
			#self.exofop_rowidx = exofop_rowidx

		try:
			#self.NEA_rowidx = int(NEA_rowidx)
			NEA_rowidx = int(NEA_rowidx)
		except:
			NEA_rowidx = np.nan 
			#self.NEA_rowidx = NEA_rowidx

		#print('final self.exofop_rowidx: ', self.exofop_rowidx)
		#print('final self.NEA_rowidx: ', self.NEA_rowidx)

		try:
			return NEA_rowidx, exofop_rowidx, NEA_targetname
		except:
			return NEA_rowidx, exofop_rowidx, row_target 




	### USER INPUT HANDLING
	elif self.telescope.lower() == 'user':
		NEA_columns = self.NEA_data.columns
		#exofop_data = pandas.read_csv(kep_fop_address, header=18)
		exofop_columns = self.exofop_data.columns
		row_known = 'y'
		#self.NEA_rowidx = np.nan 
		#self.exofop_rowidx = np.nan
		NEA_rowidx = np.nan
		exofop_rowidx = np.nan 


		try:
			return NEA_rowidx, exofop_rowidx, NEA_targetname
		except:
			return NEA_rowidx, exofop_rowidx, row_target 		














def get_properties(self, locate_neighbor='n'):
	print("calling _mp_attributes.py/get_properties()...")

	if self.telescope.lower() == 'kepler':
		if self.newlc == 'y':
			NEA_rowidx, NEA_targetname = self.find_planet_row(row_known='y')
		elif self.newlc == 'n': 
			NEA_rowidx, NEA_targetname = self.find_planet_row(row_known='n')	
		print('NEA_rowidx = ', NEA_rowidx)


	#elif self.telescope.lower() == 'k2' or self.telescope.lower() == 'tess':
	elif self.telescope.lower() == 'k2':
		if self.newlc == 'y':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='y')
		elif self.newlc == 'n':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='n')

		print('NEA_rowidx = ', NEA_rowidx)
		print('exofop_rowidx = ', exofop_rowidx)

	elif self.telescope.lower() == 'tess':
		if self.newlc == 'y':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='y')
		elif self.newlc == 'n':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='n')	

		print('[gp] NEA_rowidx = ', NEA_rowidx)
		print('[gp] exofop_rowidx = ', exofop_rowidx)


	elif self.telescope.lower() == 'user':
		target_period, target_period_uperr, target_period_lowerr = self.period, 0, 0
		target_tau0, target_tau0_uperr, target_tau0_lowerr = self.tau0, 0, 0
		target_impact, target_impact_uperr, target_impact_lowerr = self.impact, 0, 0
		target_duration, target_duration_uperr, target_duration_lowerr = self.duration_hours, 0, 0
		target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = self.rprstar, 0, 0
		target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = self.sma_AU, 0, 0
		target_insol, target_insol_uperr, target_insol_lowerr = np.nan, 0, 0
		target_rp, target_rp_uperr, target_rp_lowerr = self.rp_rearth, 0, 0


	### now with the rowidx we can access the other properties we want!
	if (self.telescope.lower() == 'kepler'):
		
		if self.target.lower().startswith('koi') or self.target.lower().startswith('kic'):
			#### THESE ARE OUTDATED
			target_period, target_period_uperr, target_period_lowerr = self.NEA_data['koi_period'][NEA_rowidx], self.NEA_data['koi_period_err1'][NEA_rowidx], self.NEA_data['koi_period_err2'][NEA_rowidx]
			target_tau0, target_tau0_uperr, target_tau0_lowerr = self.NEA_data['koi_time0bk'][NEA_rowidx], self.NEA_data['koi_time0bk_err1'][NEA_rowidx], self.NEA_data['koi_time0bk_err2'][NEA_rowidx]
			target_impact, target_impact_uperr, target_impact_lowerr = self.NEA_data['koi_impact'][NEA_rowidx], self.NEA_data['koi_impact_err1'][NEA_rowidx], self.NEA_data['koi_impact_err2'][NEA_rowidx]
			target_duration, target_duration_uperr, target_duration_lowerr = self.NEA_data['koi_duration'][NEA_rowidx], self.NEA_data['koi_duration_err1'][NEA_rowidx], self.NEA_data['koi_duration_err2'][NEA_rowidx]
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = self.NEA_data['koi_ror'][NEA_rowidx], self.NEA_data['koi_ror_err1'][NEA_rowidx], self.NEA_data['koi_ror_err2'][NEA_rowidx]
			target_rp, target_rp_uperr, target_rp_lowerr = self.NEA_data['koi_prad'][NEA_rowidx], self.NEA_data['koi_prad_err1'][NEA_rowidx], self.NEA_data['koi_prad_err2'][NEA_rowidx]
			target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = self.NEA_data['koi_sma'][NEA_rowidx], self.NEA_data['koi_sma_err1'][NEA_rowidx], self.NEA_data['koi_sma_err2'][NEA_rowidx]
			target_insol, target_insol_uperr, target_insol_lowerr = self.NEA_data['koi_insol'][NEA_rowidx], self.NEA_data['koi_insol_err1'][NEA_rowidx], self.NEA_data['koi_insol_err2'][NEA_rowidx]
			target_ldm1, target_ldm2 = self.NEA_data['koi_ldm_coeff1'][NEA_rowidx], self.NEA_data['koi_ldm_coeff2'][NEA_rowidx]
			target_eccen, target_eccen_uperr, target_eccen_lowerr = self.NEA_data['koi_eccen'][NEA_rowidx], self.NEA_data['koi_eccen_err1'][NEA_rowidx], self.NEA_data['koi_eccen_err2'][NEA_rowidx]
			target_longp, target_longp_uperr, target_longp_lowerr = self.NEA_data['koi_longp'][NEA_rowidx], self.NEA_data['koi_longp_err1'][NEA_rowidx], self.NEA_data['koi_longp_err2'][NEA_rowidx]
			target_incl, target_incl_uperr, target_incl_lowerr = self.NEA_data['koi_incl'][NEA_rowidx], self.NEA_data['koi_incl_err1'][NEA_rowidx], self.NEA_data['koi_incl_err2'][NEA_rowidx]	
			target_smass, target_smass_uperr, target_smass_lowerr = self.NEA_data['koi_smass'][NEA_rowidx], self.NEA_data['koi_smass_err1'][NEA_rowidx], self.NEA_data['koi_smass_err2'][NEA_rowidx]
			#target_snr = self.NEA_data['koi_model_snr']
			
		elif self.target.lower().startswith('kepler'):
			target_period, target_period_uperr, target_period_lowerr = self.NEA_data['pl_orbper'][NEA_rowidx], self.NEA_data['pl_orbpererr1'][NEA_rowidx], self.NEA_data['pl_orbpererr2'][NEA_rowidx]
			target_tau0, target_tau0_uperr, target_tau0_lowerr = self.NEA_data['pl_tranmid'][NEA_rowidx], self.NEA_data['pl_tranmiderr1'][NEA_rowidx], self.NEA_data['pl_tranmiderr2'][NEA_rowidx]
			target_impact, target_impact_uperr, target_impact_lowerr = self.NEA_data['pl_imppar'][NEA_rowidx], self.NEA_data['pl_impparerr1'][NEA_rowidx], self.NEA_data['pl_impparerr2'][NEA_rowidx]
			target_duration, target_duration_uperr, target_duration_lowerr = self.NEA_data['pl_trandur'][NEA_rowidx], self.NEA_data['pl_trandurerr1'][NEA_rowidx], self.NEA_data['pl_trandurerr2'][NEA_rowidx]
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = self.NEA_data['pl_ratror'][NEA_rowidx], self.NEA_data['pl_ratrorerr1'][NEA_rowidx], self.NEA_data['pl_ratrorerr2'][NEA_rowidx]
			target_rp, target_rp_uperr, target_rp_lowerr = self.NEA_data['pl_rade'][NEA_rowidx], self.NEA_data['pl_radeerr1'][NEA_rowidx], self.NEA_data['pl_radeerr2'][NEA_rowidx]
			NEA_target_a_rstar, NEA_target_a_rstar_uperr, NEA_target_a_rstar_lowerr = np.nanmedian(self.NEA_data['pl_ratror'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr2'][NEA_rowidx])
			NEA_target_Teq, NEA_target_Teq_uperr, NEA_target_Teq_lowerr = np.nanmedian(self.NEA_data['pl_eqt'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_eqterr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_eqterr2'][NEA_rowidx])
			target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = self.NEA_data['pl_orbsmax'][NEA_rowidx], self.NEA_data['pl_orbsmaxerr1'][NEA_rowidx], self.NEA_data['pl_orbsmaxerr2'][NEA_rowidx]
			target_insol, target_insol_uperr, target_insol_lowerr = self.NEA_data['pl_insol'][NEA_rowidx], self.NEA_data['pl_insolerr1'][NEA_rowidx], self.NEA_data['pl_insolerr2'][NEA_rowidx]
			#target_ldm1, target_ldm2 = np.nan, np.nan #### will need to define these elsewhere!
			#target_ldm1, target_ldm2 = self.NEA_data['koi_ldm_coeff1'][NEA_rowidx], self.NEA_data['koi_ldm_coeff2'][NEA_rowidx]
			target_eccen, target_eccen_uperr, target_eccen_lowerr = self.NEA_data['pl_orbeccen'][NEA_rowidx], self.NEA_data['pl_orbeccenerr1'][NEA_rowidx], self.NEA_data['pl_orbeccenerr2'][NEA_rowidx]
			#target_longp, target_longp_uperr, target_longp_lowerr = self.NEA_data['koi_longp'][NEA_rowidx], self.NEA_data['koi_longp_err1'][NEA_rowidx], self.NEA_data['koi_longp_err2'][NEA_rowidx]
			target_argp, target_argp_uperr, target_argp_lowerr = self.NEA_data['pl_orblper'][NEA_rowidx], self.NEA_data['pl_orblpererr1'][NEA_rowidx], self.NEA_data['pl_orbpererr2'][NEA_rowidx]
			target_incl, target_incl_uperr, target_incl_lowerr = self.NEA_data['pl_orbincl'][NEA_rowidx], self.NEA_data['pl_orbinclerr1'][NEA_rowidx], self.NEA_data['pl_orbinclerr2'][NEA_rowidx]	
			target_smass, target_smass_uperr, target_smass_lowerr = self.NEA_data['st_mass'][NEA_rowidx], self.NEA_data['st_masserr1'][NEA_rowidx], self.NEA_data['st_masserr2'][NEA_rowidx]
			target_teff, target_teff_uperr, target_teff_lowerr = self.NEA_data['st_teff'][NEA_rowidx], self.NEA_data['st_tefferr1'][NEA_rowidx], self.NEA_data['st_tefferr2'][NEA_rowidx]
			target_metal, target_metal_uperr, target_metal_lowerr = self.NEA_data['st_met'][NEA_rowidx], self.NEA_data['st_meterr1'][NEA_rowidx], self.NEA_data['st_meterr2'][NEA_rowidx]
			target_logg, target_logg_uperr, target_logg_lowerr = self.NEA_data['st_logg'][NEA_rowidx], self.NEA_data['st_loggerr1'][NEA_rowidx], self.NEA_data['st_loggerr2'][NEA_rowidx]
			target_ldm1, target_ldm2 = DKS_best_LDCmatch(Teff=target_teff, Logg=target_logg, MH=target_metal)



			#target_snr = self.NEA_data['koi_model_snr']



	elif (self.telescope.lower() == 'k2'):
		target_period, target_period_uperr, target_period_lowerr = np.nanmedian(self.NEA_data['pl_orbper'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_orbpererr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_orbpererr2'][NEA_rowidx])
		target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nanmedian(self.NEA_data['pl_tranmid'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_tranmiderr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_tranmiderr2'][NEA_rowidx])
		target_impact, target_impact_uperr, target_impact_lowerr = np.nanmedian(self.NEA_data['pl_imppar'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_impparerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_impparerr2'][NEA_rowidx])
		target_duration, target_duration_uperr, target_duration_lowerr = np.nanmedian(self.NEA_data['pl_trandur'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_trandurerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_trandurerr2'][NEA_rowidx])
		target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.nanmedian(self.NEA_data['pl_ratror'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr2'][NEA_rowidx])
		target_rp, target_rp_uperr, target_rp_lowerr = np.nanmedian(self.NEA_data['pl_rade'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_radeerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_radeerr2'][NEA_rowidx])
		target_a_rstar, target_a_rstar_uperr, target_a_rstar_lowerr = np.nanmedian(self.NEA_data['pl_ratror'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr2'][NEA_rowidx])
		target_Teq, target_Teq_uperr, target_Teq_lowerr = np.nanmedian(self.NEA_data['pl_eqt'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_eqterr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_eqterr2'][NEA_rowidx])


	elif (self.telescope.lower() == "tess"):
		if (np.isfinite(exofop_rowidx) == True) or (np.isfinite(exofop_rowidx) == np.array([True])):
			exofop_entry_present = 'y'
			print('searching for target parameters in the exofop database.')
			print(' ')
			exofop_target_period, exofop_target_period_uperr, exofop_target_period_lowerr = np.array(self.exofop_data['Period (days)'])[exofop_rowidx], np.array(self.exofop_data['Epoch (BJD) err'])[exofop_rowidx], np.array(self.exofop_data['Epoch (BJD) err'])[exofop_rowidx]
			exofop_target_impact, exofop_target_impact_uperr, exofop_target_impact_lowerr = np.nan, np.nan, np.nan 
			exofop_target_duration, exofop_target_duration_uperr, exofop_target_duration_lowerr = np.array(self.exofop_data['Duration (hours)'])[exofop_rowidx], np.array(self.exofop_data['Duration (hours) err'])[exofop_rowidx], np.array(self.exofop_data['Duration (hours) err'])[exofop_rowidx]
			exofop_target_rprstar, exofop_target_rprstar_uperr, exofop_target_rprstar_lowerr = np.sqrt(1e-6*np.array(self.exofop_data['Depth (ppm)'])[exofop_rowidx]), np.sqrt(1e-6*np.array(self.exofop_data['Depth (ppm) err'])[exofop_rowidx]), np.sqrt(1e-6*np.array(self.exofop_data['Depth (ppm) err'])[exofop_rowidx])
			exofop_target_rp, exofop_target_rp_uperr, exofop_target_rp_lowerr = np.array(self.exofop_data['Planet Radius (R_Earth)'])[exofop_rowidx], np.array(self.exofop_data['Planet Radius (R_Earth) err'])[exofop_rowidx], np.array(self.exofop_data['Planet Radius (R_Earth) err'])[exofop_rowidx]
			exofop_target_tau0, exofop_target_tau0_uperr, exofop_target_tau0_lowerr = np.array(self.exofop_data['Epoch (BJD)'])[exofop_rowidx], np.array(self.exofop_data['Epoch (BJD) err'])[exofop_rowidx], np.array(self.exofop_data['Epoch (BJD) err'])[exofop_rowidx]
		else:
			exofop_entry_present = 'n'


		if (np.isfinite(NEA_rowidx) == True) or (np.isfinite(NEA_rowidx) == np.array([True])):
			NEA_entry_present = 'y'
			NEA_target_period, NEA_target_period_uperr, NEA_target_period_lowerr = np.nanmedian(self.NEA_data['pl_orbper'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_orbpererr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_orbpererr2'][NEA_rowidx])
			NEA_target_tau0, NEA_target_tau0_uperr, NEA_target_tau0_lowerr = np.nanmedian(self.NEA_data['pl_tranmid'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_tranmiderr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_tranmiderr2'][NEA_rowidx])
			NEA_target_impact, NEA_target_impact_uperr, NEA_target_impact_lowerr = np.nanmedian(self.NEA_data['pl_imppar'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_impparerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_impparerr2'][NEA_rowidx])
			NEA_target_duration, NEA_target_duration_uperr, NEA_target_duration_lowerr = np.nanmedian(self.NEA_data['pl_trandur'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_trandurerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_trandurerr2'][NEA_rowidx])
			NEA_target_rprstar, NEA_target_rprstar_uperr, NEA_target_rprstar_lowerr = np.nanmedian(self.NEA_data['pl_ratror'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr2'][NEA_rowidx])
			NEA_target_rp, NEA_target_rp_uperr, NEA_target_rp_lowerr = np.nanmedian(self.NEA_data['pl_rade'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_radeerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_radeerr2'][NEA_rowidx])
			NEA_target_a_rstar, NEA_target_a_rstar_uperr, NEA_target_a_rstar_lowerr = np.nanmedian(self.NEA_data['pl_ratror'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_ratrorerr2'][NEA_rowidx])
			NEA_target_Teq, NEA_target_Teq_uperr, NEA_target_Teq_lowerr = np.nanmedian(self.NEA_data['pl_eqt'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_eqterr1'][NEA_rowidx]), np.nanmedian(self.NEA_data['pl_eqterr2'][NEA_rowidx])
			target_sma_AU, target_sma_uperr_AU, target_sma_lowerr_AU = self.NEA_data['pl_orbsmax'][NEA_rowidx], self.NEA_data['pl_orbsmaxerr1'][NEA_rowidx], self.NEA_data['pl_orbsmaxerr2'][NEA_rowidx]
			target_insol, target_insol_uperr, target_insol_lowerr = self.NEA_data['pl_insol'][NEA_rowidx], self.NEA_data['pl_insolerr1'][NEA_rowidx], self.NEA_data['pl_insolerr2'][NEA_rowidx]
			#target_ldm1, target_ldm2 = np.nan, np.nan #### will need to define these elsewhere!
			#target_ldm1, target_ldm2 = self.NEA_data['koi_ldm_coeff1'][NEA_rowidx], self.NEA_data['koi_ldm_coeff2'][NEA_rowidx]
			target_eccen, target_eccen_uperr, target_eccen_lowerr = self.NEA_data['pl_orbeccen'][NEA_rowidx], self.NEA_data['pl_orbeccenerr1'][NEA_rowidx], self.NEA_data['pl_orbeccenerr2'][NEA_rowidx]
			#target_longp, target_longp_uperr, target_longp_lowerr = self.NEA_data['koi_longp'][NEA_rowidx], self.NEA_data['koi_longp_err1'][NEA_rowidx], self.NEA_data['koi_longp_err2'][NEA_rowidx]
			target_argp, target_argp_uperr, target_argp_lowerr = self.NEA_data['pl_orblper'][NEA_rowidx], self.NEA_data['pl_orblpererr1'][NEA_rowidx], self.NEA_data['pl_orbpererr2'][NEA_rowidx]
			target_incl, target_incl_uperr, target_incl_lowerr = self.NEA_data['pl_orbincl'][NEA_rowidx], self.NEA_data['pl_orbinclerr1'][NEA_rowidx], self.NEA_data['pl_orbinclerr2'][NEA_rowidx]	
			target_smass, target_smass_uperr, target_smass_lowerr = self.NEA_data['st_mass'][NEA_rowidx], self.NEA_data['st_masserr1'][NEA_rowidx], self.NEA_data['st_masserr2'][NEA_rowidx]
			target_teff, target_teff_uperr, target_teff_lowerr = self.NEA_data['st_teff'][NEA_rowidx], self.NEA_data['st_tefferr1'][NEA_rowidx], self.NEA_data['st_tefferr2'][NEA_rowidx]
			target_metal, target_metal_uperr, target_metal_lowerr = self.NEA_data['st_met'][NEA_rowidx], self.NEA_data['st_meterr1'][NEA_rowidx], self.NEA_data['st_meterr2'][NEA_rowidx]
			target_logg, target_logg_uperr, target_logg_lowerr = self.NEA_data['st_logg'][NEA_rowidx], self.NEA_data['st_loggerr1'][NEA_rowidx], self.NEA_data['st_loggerr2'][NEA_rowidx]
			target_ldm1, target_ldm2 = Claret_best_LDCmatch(Teff=target_teff, Logg=target_logg, MH=target_metal)




		else:
			NEA_entry_present = 'n'


		if ((exofop_entry_present == 'y') and (NEA_entry_present == 'y')) or ((exofop_entry_present == 'n') and (NEA_entry_present == 'y')):
			#### prioritize the NEA entries
			print("USING / PRIORITIZING NEA PARAMETERS.")
			target_period, target_period_uperr, target_period_lowerr = NEA_target_period, NEA_target_period_uperr, NEA_target_period_lowerr
			target_impact, target_impact_uperr, target_impact_lowerr = NEA_target_impact, NEA_target_impact_uperr, NEA_target_impact_lowerr 
			target_duration, target_duration_uperr, target_duration_lowerr = NEA_target_duration, NEA_target_duration_uperr, NEA_target_duration_lowerr
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = NEA_target_rprstar, NEA_target_rprstar_uperr, NEA_target_rprstar_lowerr
			target_rp, target_rp_uperr, target_rp_lowerr = NEA_target_rp, NEA_target_rp_uperr, NEA_target_rp_lowerr
			target_tau0, target_tau0_uperr, target_tau0_lowerr = NEA_target_tau0, NEA_target_tau0_uperr, NEA_target_tau0_lowerr
			target_a_rstar, target_a_rstar_uperr, target_a_rstar_lowerr = NEA_target_a_rstar, NEA_target_a_rstar_uperr, NEA_target_a_rstar_lowerr
			target_Teq, target_Teq_uperr, target_Teq_lowerr = NEA_target_Teq, NEA_target_Teq_uperr, NEA_target_Teq_lowerr





		elif (exofop_entry_present == 'y') and (NEA_entry_present == 'n'):
			print("USING EXOFOP PARAMETERS.")
			target_period, target_period_uperr, target_period_lowerr = exofop_target_period, exofop_target_period_uperr, exofop_target_period_lowerr
			target_impact, target_impact_uperr, target_impact_lowerr = exofop_target_impact, exofop_target_impact_uperr, exofop_target_impact_lowerr 
			target_duration, target_duration_uperr, target_duration_lowerr = exofop_target_duration, exofop_target_duration_uperr, exofop_target_duration_lowerr
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = exofop_target_rprstar, exofop_target_rprstar_uperr, exofop_target_rprstar_lowerr
			target_rp, target_rp_uperr, target_rp_lowerr = exofop_target_rp, exofop_target_rp_uperr, exofop_target_rp_lowerr
			target_tau0, target_tau0_uperr, target_tau0_lowerr = exofop_target_tau0, exofop_target_tau0_uperr, exofop_target_tau0_lowerr
			#target_a_rstar, target_a_rstar_uperr, target_a_rstar_lowerr = exofop_target_a_rstar, exofop_target_a_rstar_uperr, exofop_target_a_rstar_lowerr
			#target_Teq, target_Teq_uperr, target_Teq_lowerr = exofop_target_Teq, exofop_target_Teq_uperr, exofop_target_Teq_lowerr		

		elif (exofop_entry_present == 'n') and (NEA_entry_present == 'n'):
			print('EXOFOP AND NEA ENTRIES BOTH MISSING. SETTING PARAMETERS TO NaN.')
			target_period, target_period_uperr, target_period_lowerr = np.nan, np.nan, np.nan
			target_impact, target_impact_uperr, target_impact_lowerr = np.nan, np.nan, np.nan
			target_duration, target_duration_uperr, target_duration_lowerr = np.nan, np.nan, np.nan
			target_rprstar, target_rprstar_uperr, target_rprstar_lowerr = np.nan, np.nan, np.nan
			target_rp, target_rp_uperr, target_rp_lowerr = np.nan, np.nan, np.nan
			target_tau0, target_tau0_uperr, target_tau0_lowerr = np.nan, np.nan, np.nan
			target_a_rstar, target_a_rstar_uperr, target_a_rstar_lowerr = np.nan, np.nan, np.nan
			target_Teq, target_Teq_uperr, target_Teq_lowerr = np.nan, np.nan, np.nan					



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


	if self.telescope.lower() == 'kepler' and (float(target_tau0) > 2454833):
		try:
			self.tau0 = float(target_tau0) - 2454833 
		except:
			self.tau0 = np.nan


	elif (self.telescope.lower() == 'k2') and (float(target_tau0) > 2454833):
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
	if self.telescope.lower() == 'kepler':
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
			#pass
			target_a_rstar = np.nan 
			self.a_rstar = np.nan 

	try:
		self.insol = float(target_insol) ### insolution in units of Earth insolation
		self.insol_err = (float(target_insol_lowerr), float(target_insol_uperr))
	except:
		#pass
		self.insol = np.nan 
		self.insol_err = np.nan 

	try:
		self.Teq = float(target_Teq) ### equilibrium temperature of the planet
		self.Teq_err = (float(target_Teq_lowerr), float(target_Teq_uperr))
	except:
		#pass
		self.Teq = np.nan 
		self.Teq_err = np.nan 


	try:
		self.eccen = float(target_eccen) ### equilibrium temperature of the planet
		self.eccen_err = (float(target_eccen_lowerr), float(target_eccen_uperr))
	except:
		#pass
		self.eccen = np.nan 
		self.eccen_err = np.nan 

	try:
		self.longp = float(target_longp) ### equilibrium temperature of the planet
		self.longp_err = (float(target_longp_lowerr), float(target_longp_uperr))
	except:
		#pass
		self.longp = np.nan 
		self.longp_err = np.nan 


	try:
		self.incl = float(target_incl) ### equilibrium temperature of the planet
		self.incl_err = (float(target_incl_lowerr), float(target_incl_uperr))
	except:
		#pass
		self.incl = np.nan 
		self.incl_err = np.nan 

	try:
		self.smass = float(target_smass)
		self.smass_err = (float(target_smass_lowerr), float(target_smass_uperr))
	except:
		#pass
		self.smass = np.nan 
		self.smass_err = np.nan 


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




def find_taus(self):
	print("calling _mp_attributes.py/find_taus().")	
	try:
		transit_midtimes = self.taus 

	except:
		try:
			transit_midtimes = [self.tau0]
			print('tau0 = ', self.tau0)
			print('looking for transit times...')
			print('self.period = ', self.period)

			if (self.period == 0.0) or np.isfinite(self.period) == False:
				manual_period_entry = input('Something wrong with the planet period... please enter a value: ')
				self.period = float(manual_period_entry)

			while (transit_midtimes[-1] - self.period) > np.nanmin(np.hstack(self.times)):
				#print('appending transit_midtime: ', transit_midtimes[-1] - self.period)
				### the transit_midtime you just added isn't the first transit!
				transit_midtimes.append(transit_midtimes[-1] - self.period)
			transit_midtimes = np.sort(transit_midtimes).tolist()
			print('appended '+str(len(transit_midtimes))+' transit midtimes.')

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
	print('calling _mp_attributes.py/mystery_solver().')
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
	print('calling _mp_attributes.py/find_neighbors().')
	print('is_neighbor = ', is_neighbor)

	if is_neighbor == 'y':
		row_known = 'n'
	else:
		row_known = 'y'
	#row_known = 'n'

	if self.telescope.lower() == 'user':
		pass 

	elif self.telescope.lower() == 'k2':
		NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known=row_known)
		if NEA_rowidx < 10:
			check_NEA_rows = np.arange(0,NEA_rowidx+11,1)
		else:
			check_NEA_rows = np.arange(NEA_rowidx-10,NEA_rowidx+10,1)

		if exofop_rowidx < 10:
			check_exofop_rows = np.arange(0,exofop_rowidx+11,1)
		else:
			check_exofop_rows = np.arange(exofop_rowidx-10,exofop_rowidx+11,1)
		print('NEA_rowidx, exofop_rowidx = ', NEA_rowidx, exofop_rowidx)
	


	elif self.telescope.lower() == 'tess':
		NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known=row_known)
		if np.isfinite(NEA_rowidx) and (NEA_rowidx < 10):
			check_NEA_rows = np.arange(0,NEA_rowidx+11,1)
		elif np.isfinite(NEA_rowidx) and (NEA_rowidx >= 10):
			check_NEA_rows = np.arange(NEA_rowidx-10,NEA_rowidx+10,1)
		else:
			check_NEA_rows = np.array([NEA_rowidx])

		if np.isfinite(exofop_rowidx) and (exofop_rowidx < 10):
			check_exofop_rows = np.arange(0,exofop_rowidx+11,1)
		
		elif np.isfinite(exofop_rowidx) and (exofop_rowidx >= 10):
			check_exofop_rows = np.arange(exofop_rowidx-10,exofop_rowidx+11,1)
		
		else:
			check_exofop_rows = np.array([exofop_rowidx])
		
		print('NEA_rowidx, exofop_rowidx = ', NEA_rowidx, exofop_rowidx)



	elif self.telescope.lower() == 'kepler':
		
		NEA_rowidx, NEA_targetname = self.find_planet_row(row_known=row_known)
		print('NEA_rowidx, NEA_targetname = ', NEA_rowidx, NEA_targetname)

		if self.target.lower().startswith('koi') or self.target.lower().startswith('kic'):
			if (type(NEA_rowidx) == list) and (len(NEA_rowidx) == 0):
					### means this object is not in the cumulative KOI list.
					print('It appears this object is not in the cumulative KOI list.')
					print('checking all '+str(len(self.NEA_data['kepoi_name']))+' MAST rows...')
					#### you're going to need to check all mast rows!
					check_NEA_rows = np.arange(0,len(self.NEA_data['kepoi_name']),1)

			elif np.isfinite(NEA_rowidx) == False:
				check_NEA_rows = np.arange(0,len(self.NEA_data['kepoi_name']),1)	

			else:
				if NEA_rowidx < 10:
					check_NEA_rows = np.arange(0,NEA_rowidx+11,1)
				else:
					check_NEA_rows = np.arange(NEA_rowidx-10,NEA_rowidx+10,1)				


		elif self.target.lower().startswith('kepler'):

			if (type(NEA_rowidx) == list) and (len(NEA_rowidx) == 0):
					### means this object is not in the cumulative KOI list.
					print('It appears this object is not in the confirmed planets list.')
					print('checking all '+str(len(self.NEA_data['pl_name']))+' MAST rows...')
					#### you're going to need to check all mast rows!
					check_NEA_rows = np.arange(0,len(self.NEA_data['pl_name']),1)

			elif np.isfinite(NEA_rowidx) == False:
				check_NEA_rows = np.arange(0,len(self.NEA_data['pl_name']),1)	


			else:
				if NEA_rowidx < 10:
					check_NEA_rows = np.arange(0,NEA_rowidx+11,1)
				else:
					check_NEA_rows = np.arange(NEA_rowidx-10,NEA_rowidx+10,1)


		else:
			if NEA_rowidx < 10:
				check_NEA_rows = np.arange(0,NEA_rowidx+11,1)
			else:
				check_NEA_rows = np.arange(NEA_rowidx-10,NEA_rowidx+10,1)

		print('NEA_rowidx = ', NEA_rowidx)


	print('checking '+str(len(check_NEA_rows))+' rows.')

	neighbor_rows = []
	neighbor_targets = []

	if self.target.lower().startswith('usr'):
		pass

	elif self.target.lower().startswith('koi'):
		print("looking for neighbors in MAST for this KOI.")

		for cr in check_NEA_rows:
			if cr <= len(self.NEA_data['kepoi_name']) - 1:
				target_term = self.NEA_targetname[:-1]
				search_term = np.array(self.NEA_data['kepoi_name'])[cr][:-1]
				if (search_term == target_term) and (cr != NEA_rowidx):
					print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.NEA_data['kepoi_name'][cr]))
					neighbor = str(self.NEA_data['kepoi_name'][cr])
					while neighbor.startswith('K') or neighbor.startswith('0'):
						neighbor = neighbor[1:]
					neighbor = 'KOI-'+str(neighbor)
					neighbor_rows.append(cr)
					neighbor_targets.append(neighbor)


	elif self.target.lower().startswith('kepler'):
		print("looking for neighbors in MAST for this Kepler target.")
		for cr in check_NEA_rows:
			#if cr <= len(self.NEA_data['kepler_name']) - 1:
			if cr <= len(self.NEA_data['pl_name']):
				target_term = self.target[7:-1]
				#search_term = np.array(self.NEA_data['kepler_name'])[cr][7:-2]
				search_term = np.array(self.NEA_data['pl_name'])[cr][7:-2]
				#if ((len(NEA_rowidx) == 0) and (search_term == target_term)) or ((search_term == target_term) and (cr != NEA_rowidx)):
				if ((search_term == target_term)) or ((search_term == target_term) and (cr != NEA_rowidx) and (search_term[-1] != NEA_target_term[-1])):				
					#print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.NEA_data['kepler_name'][cr]))
					print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.NEA_data['pl_name'][cr]))
					#neighbor = str(self.NEA_data['kepler_name'][cr])
					neighbor = str(self.NEA_data['pl_name'][cr])
					if ' ' in neighbor:
						neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
					neighbor_rows.append(cr)
					neighbor_targets.append(neighbor)

	elif self.target.lower().startswith('k2'):
		print("looking for neighbors in MAST for this K2 target.")
		for cr in check_NEA_rows:
			if cr <= len(self.NEA_data['pl_name']) - 1:
				target_term = NEA_targetname[3:-2]
				search_term = np.array(self.NEA_data['pl_name'])[cr][3:-2]
				if (search_term == target_term) and (np.array(self.NEA_data['pl_name'])[cr] != NEA_targetname):
					print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.NEA_data['pl_name'][cr]))
					neighbor = str(self.NEA_data['pl_name'][cr])
					if ' ' in neighbor:
						neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
					neighbor_rows.append(cr)
					neighbor_targets.append(neighbor)		

	else:
		if (self.telescope.lower() == "tess"):

			if self.target.lower().startswith("tic"):
				print("looking for neighbors in exofop for this TIC target.")
				for cr in check_exofop_rows:
					if cr <= len(self.exofop_data['TIC ID']) - 1:
						if (str(np.array(self.exofop_data['TIC ID'])[cr])[3:] == str(NEA_targetname)[3:]) and (cr != exofop_rowidx):
							print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.exofop_data['TIC ID'][cr]))
							neighbor = str(self.exofop_data['TIC ID'][cr])
							if ' ' in neighbor:
								neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
							neighbor_rows.append(cr)
							neighbor_targets.append('TIC '+str(neighbor))

			elif self.target.lower().startswith('toi'):
				print('looking for neighbors in exofop for this TOI.')
				for cr in check_exofop_rows:
					if cr <= len(self.exofop_data['TOI']) - 1:
						if (str(np.array(self.exofop_data['TOI'])[cr])[:-2] == str(NEA_targetname)[4:-2]) and (cr != exofop_rowidx):
							print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.exofop_data['TOI'][cr]))
							neighbor = str(self.exofop_data['TOI'][cr])
							if ' ' in neighbor:
								neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
							neighbor_rows.append(cr)
							neighbor_targets.append('TOI-'+str(neighbor))						

			else:
				try:
					print('looking fore neighbors in MAST for this TESS target.')
					for cr in check_NEA_rows:
						if cr <= len(self.NEA_data['pl_name']) - 1:
							if (np.array(self.NEA_data['pl_name'])[cr][3:-2] == NEA_targetname[3:-2]) and (np.array(self.NEA_data['pl_name'])[cr] != NEA_targetname):
								print('FOUND A NEIGHBOR FOR '+str(NEA_targetname)+': '+str(self.NEA_data['pl_name'][cr]))
								neighbor = str(self.NEA_data['pl_name'][cr])
								if ' ' in neighbor:
									neighbor = neighbor[:-2]+neighbor[-1] ### removes the space
								neighbor_rows.append(cr)
								neighbor_targets.append(neighbor)		
				except:
					pass							
	try:
		neighbor_targets.remove(self.target)
	except:
		pass

	self.neighbors = neighbor_targets









def find_TTVs(self, show_plot='n', yvar='OCmins', mask_multiple=None):
	print('calling _mp_attributes.py/find_TTVs().')
	#if mask_multiple == None:
	#	mask_multiple = self.mask_multiple 

	OC_mins = []
	OC_days = []
	OC_over_durs = []
	epochs = []

	nquarters = len(self.quarters)
	for qidx in np.arange(0,quarters,1):
		if nquarters != 1:
			qtimes, qfluxes, qerrors = self.times[qidx], self.fluxes_detrend[qidx], self.errors_detrend[qidx]
		elif nquarters == 1:
			qtimes, qfluxes, qerrors = self.times, self.fluxes_detrend, self.errors_detrend 
		for epoch, tau in enumerate(self.taus):
			if tau >= np.nanmin(qtimes) and tau <= np.nanmax(qtimes):
				### this tau is in the quarter
				### grab times up to 3x the transit duration on either side of the the tau.
				transit_idxs = np.where((qtimes >= tau - ((self.mask_multiple/2)*self.duration_days)) & (qtimes <= tau + ((self.mask_multiple/2)*self.duration_days)))[0]
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
	print('calling _mp_attributes.py/get_future_transits().')
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

	return future_taus, output_taus




