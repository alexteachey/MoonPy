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
import warnings
from datetime import datetime 
import inspect


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
from _mp_visuals import fold 
import subprocess


warnings.filterwarnings("ignore", category=UserWarning)


plt.rcParams["font.family"] = 'serif'

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/_mp_attributes.py')]


#### test the moonpy import right here
#print("TESTING THE MOONPY IMPORT (_mp_attributes.py)")
#dummy = MoonpyLC(targetID='Kepler-1625b')


def find_transit_quarters(self, times=None, fluxes=None, errors=None, locate_neighbor='n'):
	print('calling _mp_attributes.py/find_transit_quarters().')
	self.get_properties(locate_neighbor=locate_neighbor, times=times, fluxes=fluxes, errors=errors) ### calls find_planet_row() within.
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



def print_attributes(self, return_vals=False):
	five_per_line = []
	attributes = np.sort(self.__dir__())

	for nattr, attr in enumerate(attributes):
		if attr.startswith('_'):
			continue
		if len(five_per_line) < 5:
			five_per_line.append(attr)
		elif len(five_per_line) == 5:
			print(five_per_line)
			five_per_line = []

	if return_vals == True:
		return attributes 








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




def get_file_age(filepath, format='days'):
	current_time = time.time()
	if os.path.exists(filepath):
		file_creation_time = os.path.getctime(filepath)
	else:
		file_creation_time = 0
		print(filepath+' NOT FOUND. SETTING CREATION TIME = 0.')
	age_seconds = current_time - file_creation_time
	age_days = age_seconds / 86400

	if format=='days':
		return age_days
	elif format=='seconds':
		return age_seconds 


def show_function_inputs(func, return_vals=False):
	args, defaults = inspect.getfullargspec(func).args, inspect.getfullargspec(func).defaults
	print(str(func))
	print('DEFAULTS: ')
	try:
		for arg, default in zip(args, defaults):
			print(str(arg)+' = '+str(default))
	except:
		print(args, defaults)

	if return_vals == True:
		return args, defaults 



def get_databases(target_prefix):

	current_time = time.time()

	#### DEFINE FILE NAMES 
	if target_prefix.lower() in ['kepler', 'koi', 'kic']:
		mission = 'kepler'
		NEA_confirmed_address = moonpydir+'/NEA_confirmed_planets.csv'		
		NEA_candidates_address = moonpydir+'/NEA_cumkois.csv'
		exofop_address = moonpydir+'/kepler_exofop_targets.csv'
		exofop_URL = 'https://exofop.ipac.caltech.edu/kepler/download_summary_csv.php?sort=koi'

	elif target_prefix.lower() in ['k2', 'epic']:
		mission = 'k2'
		NEA_confirmed_address = moonpydir+'/NEA_confirmed_planets.csv'
		NEA_candidates_address = moonpydir+'/NEA_cumk2ois.csv'
		exofop_address = moonpydir+'/k2_exofop_targets.csv'
		exofop_URL = 'https://exofop.ipac.caltech.edu/k2/download_summary_csv.php?camp=All&sort=target'


	elif target_prefix.lower() in ['toi', 'tic', 'tess']:
		mission = 'tess'
		NEA_confirmed_address = moonpydir+'/NEA_confirmed_planets.csv'		
		NEA_candidates_address = moonpydir+'/cumtois.csv' ### doesn't exist, but whatever
		exofop_address = moonpydir+'/tess_exofop_targets.csv'		
		exofop_URL = 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv'	

	old_NEA_candidates_address = NEA_candidates_address[:-4]+'_OLD.csv'
	old_NEA_confirmed_address = NEA_confirmed_address[:-4]+'_OLD.csv'
	old_exofop_address = exofop_address[:-4]+'_OLD.csv'


	##### CHECK FILE CREATION TIMES
	NEA_confirmed_age_days = get_file_age(NEA_confirmed_address)	
	NEA_candidates_age_days = get_file_age(NEA_candidates_address)
	exofop_age_days = get_file_age(exofop_address)


	#### determine if any are out of date..
	if NEA_confirmed_age_days > 1:
		NEA_confirmed_ood = True
		print('NASA Exoplanet Archive confirmed planets database is '+str(NEA_confirmed_age_days)+' days out of date.')

	else:
		NEA_confirmed_ood = False

	if NEA_candidates_age_days > 1:
		print('NASA Exoplanet Archive candidates database is '+str(NEA_candidates_age_days)+' days out of date.')		
		NEA_candidates_ood = True
	else:
		NEA_candidates_ood = False

	if exofop_age_days > 1:
		print('ExoFOP database is '+str(exofop_age_days)+' days out of date.')			
		exofop_ood = True
	else:
		exofop_ood = False 


	opt_out_filepath = moonpydir+'/'+target_prefix+'_updated_databases_opt_out.txt'
	opt_out_age_days = get_file_age(opt_out_filepath)

	if np.any((NEA_confirmed_ood, NEA_candidates_ood, exofop_ood)): ### one day old

		#### check to see if the updated_databases_decision file exists.
		if os.path.exists(opt_out_filepath):

			#### check to see if the file is more than one day old.
			if opt_out_age_days > 1:
				#### remove it.
				os.system('rm '+opt_out_filepath)

				#### CHECK WITH USER ABOUT RE-DOWNLOADING.
				print(' ')
				print('Planet databases are missing or more than one day out of date. Do you want to download the new version? ')
				print('NOTE: this can take several minutes, and may be unnecessary if the databases are not being updated.')
				print(' ')
				download_new_csvs = input('Download updated databases? y/n: ')

			else:
				download_new_csvs = 'n'

		else:
			#### CHECK WITH USER ABOUT RE-DOWNLOADING.
			print(' ')
			print('Planet databases are missing or more than one day out of date. Do you want to download the new version? ')
			print('NOTE: this can take several minutes, and may be unnecessary if the databases are not being updated.')
			print(' ')
			download_new_csvs = input('Download updated databases? y/n: ')



	else:
		download_new_csvs = 'n'
		opt_out_age_days = get_file_age(opt_out_filepath)

	if (download_new_csvs == 'n') and (opt_out_age_days > 1): ##### don't need to be re-writing this every time you run a light curve!
		#### write this to a file so you don't have to keep asking every time.
		updated_databases_decision_file = open(moonpydir+'/'+target_prefix+'_updated_databases_opt_out.txt', mode='w')
		now = datetime.now()
		dt_string = now.strftime('%d/%m/%Y %H:%M:%S')
		updated_databases_decision_file.write('You opted not to update the databases on '+dt_string+'\n')
		updated_databases_decision_file.write('NEA confirmed planets database age is '+str(NEA_confirmed_age_days)+' days. \n')
		updated_databases_decision_file.write('NEA candidate planets database age is '+str(NEA_candidates_age_days)+' days. \n')		
		updated_databases_decision_file.write('ExoFOP planet database age is '+str(exofop_age_days)+' days. \n')		
		updated_databases_decision_file.close()



	if download_new_csvs == 'y':
		#### copy older files to make sure you have a working back-up (in case of corrupt download)
		os.system('cp '+NEA_candidates_address+' '+old_NEA_candidates_address)
		os.system('cp '+NEA_confirmed_address+' '+old_NEA_confirmed_address)	
		os.system('cp '+exofop_address+' '+old_exofop_address)

		#### access exofop credentials
		exofop_username, exofop_pw = get_exofop_credentials()


		##### load the list of desired columns
		NEA_confirmed_desired_columns = np.load(moonpydir+'/NEA_confirmed_planets_desired_columns.npy')
		print('number of desired columns: ', len(NEA_confirmed_desired_columns))

		NEA_confirmed_desired_cols = ''
		for ndc,dc in enumerate(NEA_confirmed_desired_columns):

			NEA_confirmed_desired_cols = NEA_confirmed_desired_cols+dc+','
		if NEA_confirmed_desired_cols[-1] == ',':
			NEA_confirmed_desired_cols = NEA_confirmed_desired_cols[:-1]



		##### DO THE DOWNLOADS

		#### FOR KEPLER, K2, and TESS
		print("DOWNLOADING NASA Exoplanet Archive confirmed planets file (once per day)...")
		if NEA_confirmed_ood:
			CONFIRMED_query = '"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+pscomppars+order+by+pl_name+asc&format=csv" -O "'+NEA_confirmed_address+'"'
			os.system('wget -v --tries=1 '+CONFIRMED_query)
			print(" ")


		if NEA_candidates_ood:
			if mission == 'kepler':
				print("DOWNLOADING NASA Exoplanet Archive candidate planets file...")
				#CANDIDATE_query =  '"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+cumulative+order+by+kepoi_name+asc&format=csv" -O "'+NEA_confirmed_address+'"'
				CANDIDATE_query = 'wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,koi_smass,koi_smass_err1,koi_smass_err2,koi_model_snr,ra,dec&order=kepoi_name&format=ascii" -O "'+NEA_candidates_address+'"'
				os.system('wget -v --tries=1 '+CANDIDATE_query)	
				print(' ')
				print(' ')

			elif mission == 'k2':
				print("DOWNLOADING NASA Exoplanet Archive candidate planets file...")
				#### ONLINE VIEW IS AT: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=k2pandc
				CANDIDATE_query =  '"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+k2pandc+order+by+pl_name+asc&format=csv" -O "'+NEA_candidates_address+'"'
				#CANDIDATE_query = 'wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,koi_smass,koi_smass_err1,koi_smass_err2,koi_model_snr,ra,dec&order=kepoi_name&format=ascii" -O "'+NEA_candidates_address+'"'
				os.system('wget -v --tries=1 '+CANDIDATE_query)		
				print(' ')
				print(' ')

			elif mission == 'tess':
				print("DOWNLOADING NASA Exoplanet Archive candidate planets file...")
				#### ONLINE VIEW IS AT: https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI
				CANDIDATE_query =  '"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+*+from+TOI+order+by+toi+asc&format=csv" -O "'+NEA_candidates_address+'"'
				#CANDIDATE_query = 'wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=TOI&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,koi_smass,koi_smass_err1,koi_smass_err2,koi_model_snr,ra,dec&order=kepoi_name&format=ascii" -O "'+NEA_candidates_address+'"'
				os.system('wget -v --tries=1 '+CANDIDATE_query)		
				print(' ')
				print(' ')			

		if exofop_ood:
			print('DOWNLOADING ExoFOP file (once per day)...')
			os.system('wget --tries=1 --user='+exofop_username+' --password='+exofop_pw+' "'+exofop_URL+'" -O "'+exofop_address+'"')	
			print(' ')
			print(' ')




	#### LOAD THE DATA IN
	try:
		confirmed_NEA_data = ascii.read(NEA_confirmed_address)
		confirmed_NEA_columns = confirmed_NEA_data.columns

		#### if this succeeds, you can delete the old versions.
		if os.path.exists(old_NEA_confirmed_address):
			os.system('rm '+old_NEA_confirmed_address)
	except:
		try:
			print(' ')
			print(' ')
			print('ERROR LOADING NEWLY DOWNLOADED DATABASE (NEA_confirmed_address). REVERTING TO PREVIOUS VERSION.')
			print(' ')
			print(' ')
			os.system('rm -rf '+NEA_confirmed_address)
			os.system('mv '+old_NEA_confirmed_address+' '+NEA_confirmed_address)
			confirmed_NEA_data = ascii.read(NEA_confirmed_address)		
			confirmed_NEA_columns = confirmed_NEA_data.columns
		except:
			confirmed_NEA_data = None
			confirmed_NEA_columns = np.array([])


	try:
		candidate_NEA_data = ascii.read(NEA_candidates_address)
		candidate_NEA_columns = candidate_NEA_data.columns

		#### if this succeeds, you can remove the old version.
		if os.path.exists(old_NEA_candidates_address):
			s.system('rm '+old_NEA_candidates_address)

	except:
		try:
			print('ERROR LOADING NEWLY DOWNLOADED DATABASE (NEA_candidate_address). REVERTING TO PREVIOUS VERSION.')
			os.system('rm -rf '+NEA_candidates_address)
			os.system('mv '+old_NEA_candidates_address+' '+NEA_candidates_address)
			candidate_NEA_data = ascii.read(NEA_candidates_address)				
			candidate_NEA_columns = candidate_NEA_data.columns
		except:
			candidate_NEA_data = None
			candidate_NEA_columns = np.array([])
			print('EXCEPTION ENCOUNTERED... could not load NEA_candidate_data.')


	try:
		if mission == 'kepler':
			exofop_data = pandas.read_csv(exofop_address, header=18)
		elif mission == 'k2':
			exofop_data = pandas.read_csv(exofop_address, header=10)
		elif mission == 'tess':
			exofop_data = pandas.read_csv(exofop_address)			
		exofop_columns = exofop_data.columns

		#### if this succeeds, you can remove the old version.
		if os.path.exists(old_exofop_address):
			os.system('rm '+old_exofop_address)

	except:
		try:
			print('EXCEPTION ENCOUNTERED... new ExoFOP file may be corrupted. Switching to old version.')
			#### clobber the new one, revert to the old one, and try again:
			os.system('rm '+exofop_address)
			os.system('mv '+old_exofop_address+' '+exofop_address)
			if mission == 'kepler':
				exofop_data = pandas.read_csv(exofop_address, header=18)
			elif mission == 'k2':
				exofop_data = pandas.read_csv(exofop_address, header=10)
			elif mission == 'tess':
				exofop_data = pandas.read_csv(exofop_address)		
			exofop_columns = exofop_data.columns	
		except:
			exofop_data = None
			exofop_columns = np.array([])			
			print('EXCEPTION ENCOUNTERED... could not load exofop_data.')


	#### need to figure out what to do with NEA_data -- confirmed or candidate?
	if (target_prefix == 'kepler'):
		#### if it's kepler, it's a confirmed planet, so use the confirmed planet database.
		NEA_data, NEA_columns = confirmed_NEA_data, confirmed_NEA_columns

	elif target_prefix in np.array(['koi', 'kic', 'k2', 'toi', 'tic', 'epic']):
		#### likely NOT confirmed, so use the candidate page.
		NEA_data, NEA_columns = candidate_NEA_data, candidate_NEA_columns

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
		kep_koi_address = moonpydir+'/NEA_cumkois.csv'
		kep_NEA_address = moonpydir+'/NEA_confirmed_planets.csv'

		old_kep_koi_address = kep_koi_address[:-4]+'_OLD.csv'
		old_kep_NEA_address = kep_NEA_address[:-4]+'_OLD.csv'		

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
		#k2_NEA_address = moonpydir+'/cumk2ois.txt'
		k2_NEA_address = moonpydir+'/NEA_confirmed_planets.csv'
		k2oi_address = moonpydir+'/cumko2is.csv'		

		print('row_known: ', row_known)

		if row_known == 'n':
			if str(row_target).lower().startswith('k2'):

				target_letter = str(row_target[-1])
				if ' ' in row_target:
					NEA_targetname = row_target
				else:
					NEA_targetname = str(row_target[:-1])+' '+str(target_letter)

				print('NEA_targetname: ', NEA_targetname)

				if "NEA_rowidx" in dir(self):
					NEA_rowidx = self.NEA_rowidx
				else:	
					NEA_rowidx = np.where(np.char.lower(self.NEA_data['pl_name']) == NEA_targetname.lower())[0]

				exofop_rowidx = np.nan

				print('(find_planet_row) k2 NEA_rowidx = ', NEA_rowidx)
				print('(find_planet_row) k2 exofop_rowidx = ', exofop_rowidx) 
				print("number of rows matching this description = ", len(NEA_rowidx))

				#### choose the default
				if len(NEA_rowidx) > 1:
					NEA_rowidx_idx = np.where(self.NEA_data['default_flag'][NEA_rowidx] == 1)[0]
					NEA_rowidx = NEA_rowidx[NEA_rowidx_idx]
					print('default NEA_rowidx = ', NEA_rowidx)


			elif str(row_target).lower().startswith('epic'):
				NEA_targetname = row_target[4:] ### just the numbers!
				if (NEA_targetname.startswith(' ')) or (NEA_targetname.startswith('-')):
					NEA_targetname = NEA_targetname[1:]
				

				print('NEA_targetname: ', NEA_targetname)
				try:
					if 'exofop_rowidx' in dir(self):
						exofop_rowidx = self.exofop_rowidx 
					else:
						exofop_rowidx = np.where(np.char.lower(np.array(self.exofop_data['EPIC ID']).astype(str)) == str(NEA_targetname).lower())[0]
				
				except:
					traceback.print_exc()
					print('unable to extract the exofop_rowidx')
					exofop_rowidx = np.nan 

				try:
					if "NEA_rowidx" in dir(self):
						NEA_rowidx = self.NEA_rowidx	
					else:	
						traceback.print_exc()		
						NEA_rowidx = np.where(np.char.lower(np.array(self.NEA_data['epic_name']).astype(str)) == str(NEA_targetname).lower())[0]


					#### choose the default
					if len(NEA_rowidx) > 1:
						NEA_rowidx_idx = np.where(self.NEA_data['default_flag'][NEA_rowidx] == 1)[0]
						NEA_rowidx = NEA_rowidx[NEA_rowidx_idx]
						print('default NEA_rowidx = ', NEA_rowidx)
				
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

		else:
			NEA_rowidx = self.NEA_rowidx
			exofop_rowidx = self.exofop_rowidx


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






def make_NEA_dict(self):
	#### this function will create a dictionary of all the values available in NASA Exoplanet Archive
	planet_dict = {}
	for col in self.NEA_columns:
		planet_dict[col] = self.NEA_data[col][self.NEA_rowidx]

		#### set attributes for this object using setattr(object, attrname, value)
		setattr(self, col, self.NEA_data[col][self.NEA_rowidx])


	self.NEA_dict = planet_dict





def make_vespa_starini(self, clobber='n'):

	if (clobber == 'y') or (clobber == True):
		os.system('rm -rf '+self.savepath+'/star.ini')

	#if os.path.exists(moonpydir+'/vespa'):
	#	pass
	#else:
	#	os.system('mkdir '+moonpydir+'/vespa')

	starini_filename = self.savepath+'/star.ini'
	starini_file = open(starini_filename, mode='w')
	starini_file.write('#provide spectroscopic properties if available\n')

	#### EFFECTIVE TEMPERATURE
	if np.isfinite(self.st_teff.__float__()):
		teff_input = self.st_teff
		teff_err_input = np.nanmax((np.abs(self.st_tefferr1), np.abs(self.st_tefferr2)))
		starini_file.write('Teff = '+str(teff_input)+', '+str(teff_err_input)+' #value, uncertainty\n')


	#### Metallicity feh
	if np.isfinite(self.st_met.__float__()):
		feh_input = self.st_met
		feh_err_input = np.nanmax((np.abs(self.st_meterr1), np.abs(self.st_meterr2)))
		starini_file.write('feh = '+str(feh_input)+', '+str(feh_err_input)+' #value, uncertainty\n')


	#### Logg 
	if np.isfinite(self.st_logg.__float__()):
		logg_input = self.st_logg
		logg_err_input = np.nanmax((np.abs(self.st_loggerr1), np.abs(self.st_loggerr2)))
		starini_file.write('logg = '+str(logg_input)+', '+str(logg_err_input)+' #value, uncertainty\n')	

	starini_file.write('\n')

	#### JMAG 
	if np.isfinite(self.sy_jmag.__float__()):
		jmag_input = self.sy_jmag
		jmag_err_input = np.nanmax((np.abs(self.sy_jmagerr1), np.abs(self.sy_jmagerr2)))
		starini_file.write('J = '+str(jmag_input)+', '+str(jmag_err_input)+' #value, uncertainty\n')

	#### HMAG
	if np.isfinite(self.sy_hmag.__float__()):
		hmag_input = self.sy_hmag
		hmag_err_input = np.nanmax((np.abs(self.sy_hmagerr1), np.abs(self.sy_hmagerr2)))
		starini_file.write('H = '+str(hmag_input)+', '+str(hmag_err_input)+' #value, uncertainty\n')		

	#### KMAG
	if np.isfinite(self.sy_kmag.__float__()):
		kmag_input = self.sy_kmag
		kmag_err_input = np.nanmax((np.abs(self.sy_kmagerr1), np.abs(self.sy_kmagerr2)))
		starini_file.write('K = '+str(kmag_input)+', '+str(kmag_err_input)+' #value, uncertainty\n')


	#### Kepler mag
	if np.isfinite(self.sy_kepmag.__float__()):
		kepmag_input = self.sy_kepmag
		kepmag_err_input = np.nanmax((np.abs(self.sy_kepmagerr1), np.abs(self.sy_kepmagerr2)))
		starini_file.write('Kepler = '+str(kmag_input)+', '+str(kmag_err_input)+' #value, uncertainty\n')				


	starini_file.close()

	print('created '+starini_filename)
	



def make_vespa_fppini(self, maxrad=12, secthresh=1e-4, clobber='n'):

	if (clobber == 'y') or (clobber == True):
		os.system('rm -rf '+self.savepath+'/fpp.ini')
	#if os.path.exists(moonpydir+'/vespa'):
	#	pass
	#else:
	#	os.system('mkdir '+moonpydir+'/vespa')


	#compute the cadence
	timediffs = []
	for i in np.arange(1,len(self.times[0]), 1):
		#### use the first quarter / sector, that should do well enough
		timediffs.append(self.times[0][i] - self.times[0][i-1])
	timediffs = np.array(timediffs)
	cadence_days = np.nanmedian(timediffs)

	fppini_filename = self.savepath+'/fpp.ini'
	fppini_file = open(fppini_filename, mode='w')

	fppini_file.write('name = '+str(self.target)+'\n')
	fppini_file.write('ra = '+str(self.ra)+'\n')
	fppini_file.write('dec = '+str(self.dec)+'\n')
	fppini_file.write('\n')
	fppini_file.write('period = '+str(self.period)+' #days\n')
	fppini_file.write('rprs = '+str(self.rprstar)+' #Rp/Rstar\n')
	fppini_file.write('cadence = '+str(cadence_days)+' #cadence [days]\n')
	if self.telescope.lower() == 'kepler' or self.telescope.lower() == 'k2':
		fppini_file.write('band = Kepler\n')
	fppini_file.write('photfile = '+self.savepath+'/photfile.csv #contains transit photometry\n')
	fppini_file.write('\n')
	fppini_file.write('[constraints]\n')
	fppini_file.write('maxrad = '+str(maxrad)+' # aperture radius [arcsec]\n')
	fppini_file.write('secthres = '+str(secthresh)+' # Maximum allowed depth of potential secondary eclipse\n')

	fppini_file.close()

	print('created '+fppini_filename)



def make_vespa_photfile(self, dmeth='cofiam', clobber='n'):

	if (clobber == 'y') or (clobber == True):
		os.system('rm -rf '+self.savepath+'/photfile.csv')

	try:
		print('fold times: ', self.fold_times)
	
	except:
		#### probably means it's not detrended yet
		self.detrend(dmeth=dmeth)

		try:
			fold_times, fold_fluxes, fold_errors = self.fold_times, self.fold_fluxes, self.fold_errors
		except:
			self.fold()
			fold_times, fold_fluxes, fold_errors = self.fold_times, self.fold_fluxes, self.fold_errors 

		#### fold times are in terms of phase, so to transfor into days from tmid, you need to multiply by the period
		fold_times_from_tmid = fold_times * self.period 

		photfile_idxs = np.where((fold_times_from_tmid >= -3) & (fold_times_from_tmid <= 3))[0]
		photfile_times, photfile_fluxes, photfile_errors = fold_times_from_tmid[photfile_idxs], fold_fluxes[photfile_idxs], fold_errors[photfile_idxs]

		photfile_name = self.savepath+'/photfile.csv'
		photfile = open(photfile_name, mode='w')
		for t,f,e in zip(photfile_times, photfile_fluxes, photfile_errors):
			photfile.write(str(t)+','+str(f)+','+str(e)+'\n')
		photfile.close()

		print('created '+photfile_name)



def run_vespa(self, clobber='n'):
	if (clobber == 'y') or (clobber == True):
		print("CLOBBERING VESPA FILES AND STARTING FROM SCATCH!")
		time.sleep(3)
		#### scrub the existing files 
		os.system('rm -rf '+self.savepath+'/mist*png')
		os.system('rm -rf '+self.savepath+'/mist*h5')
		os.system('rm -rf '+self.savepath+'/starfit.log')
		os.system('rm -rf '+self.savepath+'/calcfpp/log')
		os.system('rm -rf '+self.savepath+'/chains') 

	self.make_vespa_starini(clobber=clobber)
	self.make_vespa_fppini(clobber=clobber)
	self.make_vespa_photfile(clobber=clobber)


	print("ATTEMPTING TO RUN VESPA....")

	try:
		print(' ')
		print('try #1')
		#os.system('cd '+self.savepath+'; calcfpp -n 1000')
		subprocess.Popen('cd '+self.savepath+'; calcfpp -n 1000')
		print(' ')
	except:
		traceback.print_exc()
		try:
			print(' ')
			print('try #2') #### THIS IS THE ONE!!!! (BUT MAYBE WE STILL NEED TRY #1 FOR STARFIT ALL)
			print(' ')
			#os.system('cd '+self.savepath+'; starfit --all . && calcfpp -n 1000')
			subprocess.Popen('cd '+self.savepath+'; starfit --all .&& calcfpp -n 1000', shell=True)
		except:
			traceback.print_exc()
			try:
				print(' ')
				print('try #3')
				print(' ')
				#os.system('cp *h5 '+self.savepath+'; cd '+self.savepath+'; starfit --all . && calcfpp -n 1000')
				subprocess.Popen('cp *h5 '+self.savepath+'; cd '+self.savepath+'; starfit --all . && calcfpp -n 1000', shell=True)
			except:
				traceback.print_exc()
				print(' ')
				print('last except.')
				print(' ')
				#os.system('starfit --all . && cd '+self.savepath+'; calcfpp -n 1000')
				subprocess.Popen('starfit --all . && cd '+self.savepath+'; calcfpp -n 1000', shell=True)







def get_properties(self, locate_neighbor='n', times=None, fluxes=None, errors=None):
	print("calling _mp_attributes.py/get_properties()...")


	try:
		self.make_NEA_dict()
	except:
		print('An Error occurred calling make_NEA_dict().')

	if self.telescope.lower() == 'kepler':
		if self.newlc == 'y':
			NEA_rowidx, NEA_targetname = self.find_planet_row(row_known='y')
		elif self.newlc == 'n': 
			NEA_rowidx, NEA_targetname = self.find_planet_row(row_known='n')	
		print('NEA_rowidx = ', NEA_rowidx)

		self.NEA_rowidx = NEA_rowidx 


	#elif self.telescope.lower() == 'k2' or self.telescope.lower() == 'tess':
	elif self.telescope.lower() == 'k2':
		if self.newlc == 'y':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='y')
		elif self.newlc == 'n':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='n')

		self.NEA_rowidx = NEA_rowidx 
		self.exofop_rowidx = exofop_rowidx

		print('NEA_rowidx = ', NEA_rowidx)
		print('exofop_rowidx = ', exofop_rowidx)

	elif self.telescope.lower() == 'tess':
		if self.newlc == 'y':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='y')
		elif self.newlc == 'n':
			NEA_rowidx, exofop_rowidx, NEA_targetname = self.find_planet_row(row_known='n')	

		self.NEA_rowidx = NEA_rowidx
		self.exofop_rowidx = exofop_rowidx 

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
		self.find_taus(times=times, fluxes=fluxes, errors=errors)
	except:
		#traceback.print_exc()
		print("UNABLE TO CALCULATE transit times. This error may be resolved later.")





def find_taus(self, times=None, fluxes=None, errors=None):
	print("calling _mp_attributes.py/find_taus().")	
	#try:
	#	transit_midtimes = self.taus 

	#except:
	try:
		transit_midtimes = [self.tau0]
		print('self in question: ', self.target)
		print('tau0 = ', self.tau0)
		print('looking for transit times...')
		print('self.period = ', self.period)

		if type(times) != type(None):
			self.times = times 

		if type(fluxes) != type(None):
			self.fluxes = fluxes

		if type(errors) != type(None):
			self.errors = errors

		#print('self.times = ', self.times)

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
		#traceback.print_exc()
		raise Exception('an exception was raised while calling find_taus(). Consider turning on traceback if you want to know more.')



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

		if np.isfinite(NEA_rowidx):
			if NEA_rowidx < 10:
				check_NEA_rows = np.arange(0,NEA_rowidx+11,1)
			else:
				check_NEA_rows = np.arange(NEA_rowidx-10,NEA_rowidx+10,1)
		else:
			check_NEA_rows = np.array([NEA_rowidx])


		if np.isfinite(exofop_rowidx):
			if exofop_rowidx < 10:
				check_exofop_rows = np.arange(0,exofop_rowidx+11,1)
			else:
				check_exofop_rows = np.arange(exofop_rowidx-10,exofop_rowidx+11,1)
		else:
			check_exofop_rows = np.array([exofop_rowidx])


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


	self.future_taus = future_taus #### JD
	self.output_taus = output_taus #### specified format

	return future_taus, output_taus






