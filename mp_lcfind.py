from __future__ import division
import numpy as np
import astropy
from astroquery.simbad import Simbad 
from astroquery.mast import Observations
import astropy.coordinates as coord 
from astropy import units as u
import os
from astropy.io import fits as pyfits 
import time
import traceback
import socket
from urllib.request import urlretrieve
import requests 
from mp_tools import * 
import copy 

#moonpydir = os.getcwd()

moonpydir = os.path.realpath(__file__)
moonpydir = moonpydir[:moonpydir.find('/mp_lcfind.py')]
curlscript_dir = moonpydir+'/tess_sector_curlscripts'

if os.path.exists(curlscript_dir) == False:
	os.system('mkdir '+curlscript_dir)
	#### move existing curlscripts into this directory!
	os.system('mv '+moonpydir+'/*curlscript.txt '+curlscript_dir+'/')


hostname = socket.gethostname()
if ('tethys' in hostname) and ('sinica' in hostname):
	#moonpydir = '/data/tethys/Documents/Software/MoonPy'
	central_data_dir = '/data/tethys/Documents/Central_Data'
elif ('Alexs-MacBook') in hostname:
	#moonpydir = '/Users/hal9000/Documents/Software/MoonPy'
	central_data_dir = '/Users/hal9000/Documents/Central_Data'
elif 'umbriel' in hostname:
	#moonpydir = '/home/cal/ateachey/Documents/MoonPy'
	central_data_dir = '/home/cal/ateachey/Documents/Central_Data'
else:
	#moonpydir = input('Please specify the MoonPy directory (or hard-code this into moonpy.py): ')
	#central_data_dir = input("Please specify a 'central data' directory (or hard-code this into moonpy.py): ")
	### store central_data within MoonPy directory
	if os.path.exists(moonpydir+'/Central_Data'):
		pass
	else:
		os.system('mkdir '+moonpydir+'/Central_Data')
	central_data_dir = moonpydir+'/Central_Data'
#print('moonpydir = ', moonpydir)
#print('Light curves will be stored in '+central_data_dir)



#### THIS IS A MAJOR REWORKING OF THE ORIGINAL mp_lcfind.py (eventually mp_lcfind_deprecated.py).
#### THE KEY IS, The Kplr and k2plr packages are old, not available through conda, and cause all kinds of problems.
#### So this short_cadenceript will get around them entirely, doing it all under the hood, the way you want it.
#### (Note that you'll want to continue structuring the light curves as kplr does it.)
#### 


#### WILL NEED TO FOLLOW THESE INSTRUCTIONS: https://astroquery.readthedocs.io/en/latest/mast/mast.html



#####################################

# XXXXXXXXXX ALIAS FINDING XXXXXXXX #

#####################################



def Simbad_query(ra, dec, coord_format='degrees', telescope='kepler'):
	print('calling mp_lcfind.py/Simbad_query().')
	if 'h' in ra or ':' in ra:
		coord_format == 'sexagesimal'
		ra = ra.replace(' ', '') ### remove the spaces
		dec = dec.replace(' ', '') ### remove the spaces

	if coord_format=='degrees':
		nearby_objects = Simbad.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), radius='0d0m5s')
	elif coord_format=='sexagesimal':
		nearby_objects = Simbad.query_region(coord.SkyCoord(str(ra)+' '+str(dec), frame='icrs'), radius='0d0m'+str(search_radius)+'s')

	print(nearby_objects)
	best_hit = nearby_objects[0][0].decode('UTF-8')
	print("best_hit = ", best_hit)


	if (telescope.lower() == 'kepler') and (best_hit.lower().startswith('koi') == False) and (best_hit.lower().startswith('kepler') == False) and (best_hit.lower().startswith("kic") == False):
		### look for aliases
		alias_query = Simbad.query_objectids(str(best_hit))['ID']
		search_idx = 0
		while search_idx < 10: ### maximum tries
			try:
				alias = alias_query[search_idx]
				if alias.lower().startswith('koi') or alias.lower().statswith('kepler') or alias.lower().startswith("kic"):
					best_hit = alias
					break
				else:
					search_idx += 1
			except:
				break 

		print("alias = ", alias)
		best_hit = alias

		if best_hit.lower().startswith('koi'):
			object_number = float(best_hit[4:]) ### strips off "KOI-"
			### it's valuable to keep the quarters separated like this, because they should be detrended separately!
			targtype = 'koi'
			object_name = "KOI-"+str(object_number)

		elif best_hit.lower().startswith('kepler'):
			object_number = str(int(best_hit[7:]))+'b' ### strips off "Kepler-", adds b because we don't care which planet it is.
			#object_number = object_number+'b'
			targtype = 'planet'
			object_name = 'Kepler-'+str(object_number)

		elif best_hit.lower().startswith("kic"):
			object_number = int(best_hit[4:])
			targtype = 'kic'
			object_name = 'KIC'+str(object_number)

		print("object_number = ", object_number)


	elif (telscope.lower() == 'tess') and (best_hit.lower().startswith('toi') == False) and (best_hit.lower().startswith('tic') == False):
		### look for aliases if TOI or TIC number was not listed to start.
		alias_query = Simbad.query_objectids(str(best_hit))['ID']
		search_idx = 0
		while search_idx < 10: ### maximum tries
			try:
				alias = alias_query[search_idx]
				print('potential alias: ', alias)
				if alias.lower().startswith('toi') or alias.lower().startswith('tic'):
					best_hit = alias
					break
				else:
					search_idx += 1
			except:
				break 

		print("alias = ", alias)
		best_hit = alias

		if best_hit.lower().startswith('toi'):
			object_number = float(best_hit[4:]) ### strips off "KOI-"
			### it's valuable to keep the quarters separated like this, because they should be detrended separately!
			targtype = 'toi'
			object_name = "TOI-"+str(object_number)


		elif best_hit.lower().startswith("tic"):
			object_number = int(best_hit[4:])
			targtype = 'tic'
			object_name = 'TIC'+str(object_number)

		try:
			print("object_number = ", object_number)
		except:
			pass


	return object_name, object_number 







def find_KIC_alias(target_name):
	print('calling mp_lcfind.py/find_KIC_alias().')
	#### FOR USE WITH KEPLER LIGHT CURVES.
	### in order to return a KIC number, you need to remove final letters and decimals.
	star_number = target_name
	if '.' in star_number:
		star_number = star_number[:star_number.find('.')] ### this will remove the decimal and final numbers
	if star_number[-1] in ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']:
		star_number = star_number[:-1]

	target_aliases = []
	alias_search_results = Simbad.query_objectids(star_number)
	for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
		target_alias = alias_search_results[alidx][0]
		target_aliases.append(target_alias)

		if target_alias.lower().startswith('kic'):
			kic_number = target_alias
			break
	try:
		return kic_number
	except:		
		return target_aliases


def find_TIC_alias(target_name):
	print('calling mp_lcfind.py/find_TIC_alias().')
	#### FOR USE WITH TESS LIGHT CURVES
	### in order to return a KIC number, you need to remove final letters and decimals.
	star_number = target_name
	if '.' in star_number:
		star_number = star_number[:star_number.find('.')] ### this will remove the decimal and final numbers
	if star_number[-1] in ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']:
		star_number = star_number[:-1]

	target_aliases = []
	try:
		alias_search_results = Simbad.query_objectids(star_number)
		print('len(alias_search_results): ', len(alias_search_results))
		for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
			target_alias = alias_search_results[alidx][0]
			target_aliases.append(target_alias)

			if target_alias.lower().startswith('tic'):
				tic_number = target_alias
				self.ticnum = tic_number
				break

	except:
		print('SIMBAD WAS UNABLE TO FIND THE TARGET.')
		try:
			target_aliases.np.array([self.ticnum])
			tic_number = self.ticnum
			print("TIC# LOCATED WITHIN EXOFOP.")
		except:
			try:
				find_planet_row()
				target_aliases.np.array([self.ticnum])
				tic_number = self.ticnum
				print("TIC# LOCATED WITHIN EXOFOP.")
			except:
				try:
					tic_number = 'TIC '+np.array([self.exofop_data[self.exofop_rowidx]])[0]
				except:
					print('TARGET ALIAS IS JUST ORIGINAL TARGET NAME INPUT.')
					target_aliases = np.array([target_name])

	try:
		return tic_number
	except:		
		return target_aliases



def find_EPIC_alias(target_name):
	print('calling mp_lcfind.py/find_EPIC_alias().')
	#### FOR USE WITH k2 LIGHT CURVES.
	### in order to return a EPIC number, you need to remove final letters and decimals.
	star_number = target_name
	if '.' in star_number:
		star_number = star_number[:star_number.find('.')] ### this will remove the decimal and final numbers
	if star_number[-1] in ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l']:
		star_number = star_number[:-1]

	target_aliases = []
	alias_search_results = Simbad.query_objectids(star_number)
	for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
		target_alias = alias_search_results[alidx][0]
		target_aliases.append(target_alias)

		if target_alias.lower().startswith('epic'):
			kic_number = target_alias
			break
	try:
		return kic_number
	except:		
		return target_aliases





#####################################

# XXXXXXXXXX URL GENERATION XXXXXXXX #

#####################################



def kepler_URL_generator(KIC, short_cadence=False):
	print('calling mp_lcfind.py/kepler_URL_generator()')
	### this function takes in a KIC number (with or without a prefix) and produces a URL and wget script. 

	"""
	Kepler data may be retrieved without submitting a batch request. 
	Both lightcurve and target pixel files are online and stored in two directory trees. 
	For the browser paths below KKKKKKKKK is the KIC ID and XXXX is the first 4 digits 
	of the KIC ID including the initial zeros.

	Lightcurves

	http://archive.stsci.edu/pub/kepler/lightcurves/XXXX/KKKKKKKKK
	http://archive.stsci.edu/pub/kepler/lightcurves/0014/001429092/
	"""

	query_format_number = KIC

	#### you need to take the kic number, remove kic, and add zeroes until length = 9.
	if query_format_number.lower().startswith('kic'):
		query_format_number = query_format_number[3:]
	while (query_format_number.startswith('-')) or (query_format_number.startswith(' ')):
		query_format_number = query_format_number[1:]
	### now it should be just a bunch of numbers:
	while len(query_format_number) < 9:
		query_format_number = '0'+query_format_number

	assert len(query_format_number) == 9

	first_four_numbers = query_format_number[:4]


	#http://archive.stsci.edu/pub/kepler/lightcurves/0047/004760478/kplr004760478-2010078095331_llc.fits
	final_URL = 'http://archive.stsci.edu/pub/kepler/lightcurves/'+first_four_numbers+'/'+query_format_number+'/'
	if os.path.exists(central_data_dir+'/Kepler_lightcurves') == False:
		os.system('mkdir '+central_data_dir+'/Kepler_lightcurves')

	download_directory = central_data_dir+'/Kepler_lightcurves/KIC'+str(query_format_number)
	if os.path.exists(download_directory):
		pass
	else:
		os.system('mkdir '+download_directory)

	if short_cadence == True:
		wget_command = "wget -q -nH --cut-dirs=6 -r -l0 -c -N -np -R --reject tar --accept fits 'index*' -erobots=off "+final_URL+" -P "+download_directory+"/"
	else:
		wget_command = "wget -q -nH --cut-dirs=6 -r -l0 -c -N -np -R --reject tar --accept llc.fits 'index*' -erobots=off "+final_URL+" -P "+download_directory+"/"



	return final_URL, wget_command, download_directory



def k2_URL_generator(EPIC):
	print('calling mp_lcfind.py/k2_URL_generator()')
	### this function takes in a EPIC number (with or without a prefix) and produces a URL and wget script. 

	"""
	https://archive.stsci.edu/k2/download_options.html

	Lightcurves

	Individual files may be downloaded via FTP or through your browser (HTTPS) to download K2 data and catalogs. 
	For FTP, connect to archive.stsci.edu anonymously and cd to pub/k2 You will see the available directories 
	using ls. For HTTP, just go to https://archive.stsci.edu/pub/k2/. Examples for the browser paths to light 
	curves and target pixel files are shown below, where XXXXYYZZZ is the EPIC ID and N is the campaign number. 

	Files are available online in subdirectories in the form:
	https://archive.stsci.edu/pub/k2/lightcurves/cN/XXXX00000/YY000. For example:
	https://archive.stsci.edu/pub/k2/lightcurves/c3/212200000/35000/

	### THEN THE FILE IS https://archive.stsci.edu/pub/k2/lightcurves/c3/212200000/35000/ktwo212235321-c03_llc.fits


	"""

	query_format_number = EPIC

	#### you need to take the kic number, remove kic, and add zeroes until length = 9.
	if query_format_number.lower().startswith('epic'):
		query_format_number = query_format_number[4:] ### take off the epic
	while (query_format_number.startswith('-')) or (query_format_number.startswith(' ')):
		query_format_number = query_format_number[1:] #### remove initial "-" or " "
	
	### now it should be just a bunch of numbers:
	#### THE ACTION BELOW IS UNNECCESSARY FOR K2, I THINK.
	#while len(query_format_number) < 9:
	#	query_format_number = '0'+query_format_number

	assert len(query_format_number) == 9

	XXXX = str(query_format_number[:4])
	YY = str(query_format_number[4:6]) 


	#final_URL = 'http://archive.stsci.edu/pub/kepler/lightcurves/'+first_four_numbers+'/'+query_format_number+'/'
	#final_URL = 'https://archive.stsci.edu/pub/k2/lightcurves/c3/212200000/35000/'

	download_directory = central_data_dir+'/K2_lightcurves/EPIC'+str(query_format_number)
	if os.path.exists(download_directory):
		pass
	else:
		os.system('mkdir '+download_directory)


	#### since you don't know which campaign it's in, generate a list of final_URLs, a list of wget_commands,
	final_URLs = []
	wget_commands = []
	filepaths = []

	for campaign_number in np.arange(0,19,1):
		print('campaign # ', campaign_number)
		if campaign_number < 10:
			filename = 'ktwo'+str(query_format_number)+'-c0'+str(campaign_number)+'_llc.fits'
			final_URL = 'https://archive.stsci.edu/pub/k2/lightcurves/c'+str(campaign_number)+'/'+XXXX+'00000/'+YY+'000/'+filename #ktwo'+str(query_format_number)+'-c0'+str(campaign_number)+'_llc.fits'
		else:
			filename = 'ktwo'+str(query_format_number)+'-c'+str(campaign_number)+'_llc.fits'
			final_URL = 'https://archive.stsci.edu/pub/k2/lightcurves/c'+str(campaign_number)+'/'+XXXX+'00000/'+YY+'000/'+filename #ktwo'+str(query_format_number)+'-c'+str(campaign_number)+'_llc.fits'
		final_URLs.append(final_URL)

		filepath = download_directory+'/'+filename
		filepaths.append(filepath)

		#wget_command = "wget -q -nH --cut-dirs=6 -r -l0 -c -N -np -R 'index*' -erobots=off "+final_URL+" -P "+download_directory+"/"
		#wget_command = 'wget -nH --cut-dirs=6 -r -l0 -c -N -np -R "index*" -erobots=off '+final_URL+' -P '+download_directory+'/'
		#os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,ra,dec&order=kepoi_name&format=ascii" -O "'+kep_mast_address+'"')
		wget_command = 'wget --tries=1 -N "'+final_URL+'" -P '+download_directory+'/'
		wget_commands.append(wget_command)


	return final_URLs, wget_commands, download_directory, filepaths










#####################################

# XXXXXXXXXX DOWNLOADING XXXXXXXX #

#####################################




def kepler_fits_download(target_name, clobber='n'):
	print('calling mp_lcfind.py/kepler_fits_downloader().')
	### first extract the KIC number
	try:
		KIC_name = find_KIC_alias(target_name)
		print('KIC alias = ', KIC_name)

		KIC_URL, KIC_wget, KIC_download_dir = functimer(kepler_URL_generator(KIC_name))
		#print("KIC wget = ", KIC_wget)

		#### check if KIC_download_dir files already exist... if they do, skip it!
		if (len(os.listdir(KIC_download_dir)) > 0) and (clobber == 'n'):
			print('light curve files already exist, no need to re-download.')
		else:
			print('wgetting KIC light curves...')
			os.system(KIC_wget)
			print('done.')

	except:
		traceback.print_exc()
		raise Exception("something went wrong in 'kepler_fits_download()' (see traceback).")


def k2_fits_download(target_name, clobber='n'):
	print('calling mp_lcfind.py/k2_fits_download().')
	try:
		EPIC_name = find_EPIC_alias(target_name)

		EPIC_URLs, EPIC_wgets, EPIC_download_dir, EPIC_filepaths = k2_URL_generator(EPIC_name) ### EPIC_URLs and EPIC_wgets are LISTS!!!!

		if (len(os.listdir(EPIC_download_dir)) > 0) and (clobber == 'n'):
			print("light curve files already exist, no need to download again.")
		
		else:
			##### build a list of the EPIC_URLs that have been tried before.
			unobserved_campaigns = []
			
			if os.path.exists(EPIC_download_dir+'/unobserved_campaigns.txt'):
				#### open it
				unobserved_campaigns_file = open(EPIC_download_dir+'/unobserved_campaigns.txt', mode='r')
				for nline,line in enumerate(unobserved_campaigns_file):
					unobserved_campaigns.append(line[:line.find('\n')]) ### make sure you have just the EPIC_URL, not those final \n's. 
				unobserved_campaigns_file.close()
			else:
				#### start a new list.
				unobserved_campaigns = []


			print('wgetting EPIC light curves...')
			for EPIC_URL, EPIC_wget, EPIC_filepath in zip(EPIC_URLs, EPIC_wgets, EPIC_filepaths):
				if EPIC_URL in unobserved_campaigns:
					print('target is known to not be observed in this campaign.')
					
				else:
					print('looking for ', EPIC_URL)
					os.system(EPIC_wget)
					print('done.')

					print('EPIC filepath: ', EPIC_filepath)
					if os.path.exists(EPIC_filepath) == False:
						unobserved_campaigns_file = open(EPIC_download_dir+'/unobserved_campaigns.txt', mode='a')
						unobserved_campaigns_file.write(EPIC_URL+'\n')
						unobserved_campaigns_file.close()
						print('added to the unobserved_campaigns.txt file.')



		#### check to see if the file exists -- if it doesn't... add it to the unobserved_campaigns file.


	except:
		traceback.print_exc()
		raise Exception("something went wrong in 'k2_fits_download()' (see traceback).")



#####################################

# XXXXXXXXXX UNPACK FITS XXXXXXXX #

#####################################



def kepler_unpack_fits(target_name, short_cadence=False, long_cadence=True):
	print('calling mp_lcfind.py/kepler_unpack_fits().')
	#### first thing you have to do is make sure the damn thing exists!
	KIC_directory = functimer(kepler_URL_generator(find_KIC_alias(target_name))[2])
	try:
		#### use kepler_fits_download() to grab the KIC_download_dir and unpack these fits.
		#KIC_directory = kepler_fits_download(target_name, download=download)
		KIC_directory_files = os.listdir(KIC_directory)
		print("KIC_directory_files = '", KIC_directory_files)

		### start a dictionary with quarters, times, fluxes, errors, flags
		kic_quarters_dict = {}

		KIC_fits_files = []

		### if short_cadence=True, you need to check whether there is a slc version!
		### if short_cadence=False, reject slcs!
		reject_wrong_cadence_files = []
		for i,ki in enumerate(KIC_directory_files):
			if ('slc.fits' in ki) and (short_cadence == False):
				reject_wrong_cadence_files.append(ki)
			if ('llc.fits' in ki) and (long_cadence == False):
				reject_wrong_cadence_files.append(ki)
			elif (short_cadence == True) and (long_cadence == True):
				##### now we need to give preference for short_cadence over long_cadence
				if 'llc.fits' in ki:
					filename_before_suffix = ki[:ki.find('_llc.fits')]
					continue

				elif 'slc.fits' in ki:
					filename_before_suffix = ki[:ki.find('_slc.fits')]
					reject_wrong_cadence_files.append(filename_before_suffix+'_llc.fits')


		print('reject wrong cadence files: ', reject_wrong_cadence_files)

		for Kdf in KIC_directory_files:
			if Kdf in reject_wrong_cadence_files:
				continue

			if ('.fits' in Kdf) and (Kdf not in reject_wrong_cadence_files):

				print('reading ', Kdf)

				KIC_fits_files.append(Kdf)

				kdf_path = KIC_directory+'/'+Kdf	
				kdf_file = pyfits.open(kdf_path)	

				header = kdf_file[0].header
				data = kdf_file[1].data
				quarter = header['QUARTER']

				kic_times = data['TIME']
				kic_sap_flux = data['SAP_FLUX']
				kic_sap_err = data['SAP_FLUX_ERR']
				kic_pdc_flux = data['PDCSAP_FLUX']
				kic_pdc_err = data['PDCSAP_FLUX_ERR']
				kic_sap_qual = data['SAP_QUALITY'] 

				#### create a dictionary JUST FOR THIS QUARTER
				kic_quarter_dict = {}
				kic_quarter_dict['TIME'] = kic_times
				kic_quarter_dict['SAP_FLUX'] = kic_sap_flux
				kic_quarter_dict['SAP_FLUX_ERR'] = kic_sap_err
				kic_quarter_dict['PDCSAP_FLUX'] = kic_pdc_flux
				kic_quarter_dict['PDCSAP_FLUX_ERR'] = kic_pdc_err
				kic_quarter_dict['SAP_QUALITY'] = kic_sap_qual 

				### now add this dictionary to the kic_quarter_dict 
				kic_quarters_dict[quarter] = kic_quarter_dict 

		### output the kic_quarters_dict 
		return kic_quarters_dict 
	except:
		traceback.print_exc()
		raise Exception("Something went wrong in 'kepler_unpack_fits()' (see traceback).")




def k2_unpack_fits(target_name):
	print('calling mp_lcfind.py/k2_unpack_fits()')

	#### first thing you have to do is make sure the damn thing exists!
	EPIC_directory = k2_URL_generator(find_EPIC_alias(target_name))[2] 
	try:
		#### use kepler_fits_download() to grab the KIC_download_dir and unpack these fits.
		#KIC_directory = kepler_fits_download(target_name, download=download)
		EPIC_directory_files = os.listdir(EPIC_directory)

		### start a dictionary with quarters, times, fluxes, errors, flags
		epic_quarters_dict = {}

		epic_fits_files = []
		for Edf in EPIC_directory_files:
			if '.fits' in Edf:
				epic_fits_files.append(Edf)

				edf_path = EPIC_directory+'/'+Edf	
				edf_file = pyfits.open(edf_path)	

				header = edf_file[0].header
				data = edf_file[1].data
				quarter = header['CAMPAIGN']

				epic_times = data['TIME']
				epic_sap_flux = data['SAP_FLUX']
				epic_sap_err = data['SAP_FLUX_ERR']
				epic_pdc_flux = data['PDCSAP_FLUX']
				epic_pdc_err = data['PDCSAP_FLUX_ERR']
				epic_sap_qual = data['SAP_QUALITY'] 

				#### create a dictionary JUST FOR THIS QUARTER
				epic_quarter_dict = {}
				epic_quarter_dict['TIME'] = epic_times
				epic_quarter_dict['SAP_FLUX'] = epic_sap_flux
				epic_quarter_dict['SAP_FLUX_ERR'] = epic_sap_err
				epic_quarter_dict['PDCSAP_FLUX'] = epic_pdc_flux
				epic_quarter_dict['PDCSAP_FLUX_ERR'] = epic_pdc_err
				epic_quarter_dict['SAP_QUALITY'] = epic_sap_qual 

				### now add this dictionary to the kic_quarter_dict 
				epic_quarters_dict[quarter] = epic_quarter_dict 

		### output the kic_quarters_dict 
		return epic_quarters_dict 
	except:
		traceback.print_exc()
		raise Exception("Something went wrong in 'k2_unpack_fits()' (see traceback).")





def kplr_target_download(targID, targtype='koi', quarters='all', lc_format='pdc', telescope='kepler', clobber='n', short_cadence=False):
	print('calling mp_lcfind.py/kplr_target_download().')
	#### using the functions developed above
	### first, try unpacking without downloading.
		
	### FOR MOONPY COMPATIBILITY, THIS REQUIRES SUPPORT FOR BOTH KEPLER AND K2 LIGHT CURVES.

	if telescope.lower() == 'kepler':
		kepler_fits_download(targID, clobber=clobber)
		try:
			kepler_lc_dictionary = kepler_unpack_fits(targID, short_cadence=short_cadence)
		
		except:
			print('first except triggered.')
			time.sleep(2)
			try:
				print('lc may not have been downloaded. Attempting to download...')
				kepler_fits_download(targID, clobber=clobber)
				### after it's been downloaded, you can try to unpack again.
				kepler_lc_dictionary = kepler_unpack_fits(targID, short_cadence=short_cadence)
			
			except:
				print('second except triggered.')
				time.sleep(2)
				traceback.print_exc()
				raise Exception("Something went wrong in 'kplr_target_download()' (see traceback).")

	elif telescope.lower() == 'k2':
		k2_fits_download(targID, clobber=clobber)
		try:
			kepler_lc_dictionary = k2_unpack_fits(targID)
		
		except:
			print('first except triggered.')
			time.sleep(2)
			
			try:
				print('lc may not have been downloaded. Attempting to download...')
				k2_fits_download(targID, clobber=clobber)
				### after it's been downloaded, you can try to unpack again.
				kepler_lc_dictionary = k2_unpack_fits(targID)
			
			except:
				print('second except triggered.')
				time.sleep(2)
				traceback.print_exc()
				raise Exception("Something went wrong in 'kplr_target_download()' (see traceback).")


	### now you should have kepler_lc_dictionary, that's ready to return to spit out in the original format.
	### first, sort the quarters

	sorted_quarters = np.sort(list(kepler_lc_dictionary.keys()))
	if quarters == 'all':
		output_quarters = np.array(sorted_quarters)
	else:
		output_quarters = quarters


	### make a bunch of lists -- this kinda sucks, but it's backwards compatible with the old kplr output.
	kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, kobj_sap_fluxes, kobj_sap_errors = [], [], [], [], [], [], []

	for oq in output_quarters:
		quarter_entry = kepler_lc_dictionary[int(oq)]

		#### need to find finite_idxs
		times_finite = np.where(np.isfinite(quarter_entry['TIME']))[0]
		pdc_flux_finite = np.where(np.isfinite(quarter_entry['PDCSAP_FLUX']))[0]
		pdc_fluxerr_finite = np.where(np.isfinite(quarter_entry['PDCSAP_FLUX_ERR']))[0]
		pdc_flags_finite = np.where(np.isfinite(quarter_entry['SAP_QUALITY']))[0]
		sap_flux_finite = np.where(np.isfinite(quarter_entry['SAP_FLUX']))[0]
		sap_fluxerr_finite = np.where(np.isfinite(quarter_entry['SAP_FLUX_ERR']))[0]


		finite_idxs, finite_idx_counts = np.unique(np.concatenate((times_finite, pdc_flux_finite, pdc_fluxerr_finite, pdc_flags_finite, sap_flux_finite, sap_fluxerr_finite)), return_counts=True)
		#### there are six arrays being concatenated... therefore, you only want the indices finite_idx_counts == 6.
		finite_idxs = finite_idxs[finite_idx_counts == 6]

		#print('len(finite_idxs) = ', len(finite_idxs))

		#raise Exception("Look what you've done!")

		kobj_times.append(quarter_entry['TIME'][finite_idxs])
		kobj_pdc_fluxes.append(quarter_entry['PDCSAP_FLUX'][finite_idxs])
		kobj_pdc_errors.append(quarter_entry['PDCSAP_FLUX_ERR'][finite_idxs])
		kobj_flags.append(quarter_entry['SAP_QUALITY'][finite_idxs])
		kobj_sap_fluxes.append(quarter_entry['SAP_FLUX'][finite_idxs])
		kobj_sap_errors.append(quarter_entry['SAP_FLUX_ERR'][finite_idxs])
		kobj_quarters.append(oq)

	kobj_times = np.array(kobj_times, dtype=object)
	kobj_pdc_fluxes = np.array(kobj_pdc_fluxes, dtype=object)
	kobj_pdc_errors = np.array(kobj_pdc_errors, dtype=object)
	kobj_flags = np.array(kobj_flags, dtype=object)
	kobj_sap_fluxes = np.array(kobj_sap_fluxes, dtype=object)
	kobj_sap_errors = np.array(kobj_sap_errors, dtype=object)
	kobj_quarters = np.array(kobj_quarters, dtype=object)



	print('kobj_quarters = ', kobj_quarters)
	print('len(kobj_quarters) = ', len(kobj_quarters))


	### it's valuable to keep the quarters separated like this, because they should be detrended separately!
	if lc_format == 'pdc':
		return kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters
	elif lc_format == 'sap':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters 
	elif lc_format == 'both':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters 





def kplr_coord_download(ra, dec, coord_format='degrees', quarters='all', search_radius=5, lc_format='pdc', clobber='n', short_cadence=False):
	print('calling mp_lcfind.py/kplr_coord_download().')
	### find the object in Simbad using it's coordinates, and call kplr_target_download

	print('ra,dec = ', ra, dec)	
	object_name, object_number = Simbad_query(ra=ra, dec=dec, coord_format=coord_format, telescope=self.telescope)

	if lc_format == 'pdc':
		kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, short_cadence=short_cadence)
	elif lc_format == 'sap':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, short_cadence=short_cadence)
	elif lc_format == 'both':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, short_cadence=short_cadence, clobber=clobber)


	### it's valuable to keep the quarters separated like this, because they should be detrended separately!
	if lc_format == 'pdc':
		return kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'sap':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'both':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name




def tess_coord_download(ra, dec, coord_format='degrees', quarters='all', search_radius=5, lc_format='pdc', short_cadence=False):
	print('calling mp_lcfind.py/tess_coord_download().')
	### find the object in Simbad using it's coordinates, and call kplr_target_download

	print('ra,dec = ', ra, dec)
	try:
		object_name, object_number = Simbad_query(ra=ra, dec=dec, coord_format=coord_format, telescope=self.telescope)
	except:
		print('SIMBAD QUERY FAILED.')
		#### try to find TIC number another way
		try:
			ticnumber = self.ticnum
			object_name, object_number = ticnumber, ticnumber[4:]
		except:
			try:
				ticnumber = np.array(self.exofop_data['TIC ID'][self.exofop_rowidx])[0]
				object_name, object_number = ticnumber, ticnumber[4:]
			except:
				print('EXHAUSTED ALL TRIES TO FIND TIC.')


	if lc_format == 'pdc':
		kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = tess_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, short_cadence=short_cadence)
	elif lc_format == 'sap':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters = tess_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, short_cadence=short_cadence)
	elif lc_format == 'both':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = tess_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, short_cadence=short_cadence)


	### it's valuable to keep the quarters separated like this, because they should be detrended separately!
	if lc_format == 'pdc':
		return kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'sap':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'both':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name



def tess_target_download(targID, sectors='all', short_cadence=True, lc_format='pdc', delete_fits='n', is_neighbor='n'):
	print('calling mp_lcfind.py/tess_target_download().')
	### this function interfaces with MASS to download light curves based on the TIC #.

	print('targID: ', targID)
	print('type(targID): ', type(targID))
	if type(targID) != str:
		targID = targID[0]

	if os.path.exists(moonpydir+'/TESS_lcs'):
		pass
	else:
		os.system('mkdir '+moonpydir+'/TESS_lcs')


	### first try a simple curl script!

	all_times = []
	all_fluxes = []
	all_errors = []
	all_flags = []
	sectors = []
	lcfiles = []

	try:
		if targID.lower().startswith('tic'):
			ticnum = str(targID)[3:]
			if ticnum.startswith(' '):
				ticnum = str(ticnum)[1:]
			elif ticnum.startswith('-'):
				ticnum = str(ticnum)[1:]
		
		else:
			#if targID.lower().startswith('toi'):
			#### find it in exofop!
			ticnum = self.exofop_data['TIC ID'][self.exofop_rowidx] #### I'm not sure where the self comes from, but it works! (for known TOIs)
			#ticnum = exofop_data['TIC ID'][exofop_rowidx]
			#ticnum = str(targID) ### you've already handled the TOI already!
			if ticnum.startswith(' '):
				ticnum = ticnum[1:]



		#### now we have the ticnum, we can establish the download directory. 
		download_directory = central_data_dir+'/TESS_lightcurves/TIC'+str(ticnum)
		if os.path.exists(download_directory):
			pass
		else:
			os.system('mkdir '+download_directory)	


		### prepare for the query
		query_num = ticnum
		while len(query_num) < 16:
			query_num = '0'+query_num
		assert len(query_num) == 16

		sector_prefixes, sector_suffixes = {}, {}

		### these can be found at archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html,
		### in the tesscurl_sector_NN_lc.sh files.

		nsectors = 99
		nactual_sectors = 0
		last_existing_sector = 1
		nempty_in_a_row = 0
		sector_numbers = []
		empty_curlscript_sectors = []


		sector_filename = 'tess_sectors.txt'
		if os.path.exists(moonpydir+'/'+sector_filename):
			#### see if it is more than one day old:
			current_time = time.time()
			sectorfile_createtime = os.path.getctime(moonpydir+'/'+sector_filename)

			if current_time - sectorfile_createtime > 86400: #### more than one day
				### we're going to make a new one
				create_new_sectorfile = True
			else:
				create_new_sectorfile = False 

		else:
			#### just create it
			create_new_sectorfile = True 
			os.system('touch '+moonpydir+'/'+sector_filename)



		if create_new_sectorfile == True:

			sectorfile = open(moonpydir+'/'+sector_filename, mode='w')

			for sector in np.arange(1,nsectors,1):
				### get the curl script... then extract the prefixes and suffixes from the first line.
				print('looking for sector: ', sector)

				try:
					if os.path.exists(curlscript_dir+'/sector'+str(sector)+"_curlscript.txt"):
						print('curlscript exists.')
						pass


					else:
						print('attempting to download curlscript for the first time.')
						sector_curl_URL = 'http://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_'+str(sector)+'_lc.sh'
						print('sector_curl_URL')
						os.system('wget --tries=1 -N "'+sector_curl_URL+'" -O '+curlscript_dir+'/sector'+str(sector)+"_curlscript.txt")

					curltxt = open(curlscript_dir+'/sector'+str(sector)+'_curlscript.txt', mode='r')
					first_line = curltxt.readline()
					second_line = curltxt.readline()
					sector_prefix = second_line[16:40]
					sector_suffix = second_line[56:71] 
					### now read the first line of that 
					sector_prefixes[sector], sector_suffixes[sector] = sector_prefix, sector_suffix
					print('sector_prefix = ', sector_prefix)
					print('len(sector_prefix) = ', len(sector_prefix))



					if len(sector_prefix) > 0:
						#print("sector_prefix, sector_suffix = ", sector_prefix, sector_suffix)
						nactual_sectors += 1
						#last_existing_sector = copy.copy(sector) 
						nsectors = sector
						### reset nempty_in_a_row
						nempty_in_a_row = 0 

						sectorfile.write('sector '+str(sector)+'\n')


					else:
						print('sector '+str(sector)+' curlscript file was empty...')
						print("REMOVING!")
						empty_curlscript_sectors.append(sector)
						nempty_in_a_row += 1
						#### remove it!
						os.system('rm -rf '+curlscript_dir+'/sector'+str(sector)+'_curlscript.txt')
						#break
						#continue
						if nempty_in_a_row == 5:
							break
						else:
							continue

				except:
					traceback.print_exc()
					break_loop = input('Do you want to break loop? y/n: ')
					if break_loop == 'y':
						break

			sectorfile.close()
			


		elif create_new_sectorfile == False:
			##### means you already queried the database once today, no need to do it again!

			sector_numbers = []

			nactual_sectors = 0 
			#nsectors = 0 
			#### just open it, so you can count the number of sectors
			sectorfile = open(moonpydir+'/'+sector_filename, mode='r')
			for nline, line in enumerate(sectorfile):
				sector_numbers.append(int(line.split()[1]))
				nactual_sectors += 1 	
				nsectors += 1		
			sectorfile.close()

			#nsectors = nactual_sectors
			#print('nsectors = ', nsectors)			

			for sector in sector_numbers:
				#print('sector '+str(sector))

				curltxt = open(curlscript_dir+'/sector'+str(sector)+'_curlscript.txt', mode='r')
				first_line = curltxt.readline()
				second_line = curltxt.readline()
				sector_prefix = second_line[16:40]
				sector_suffix = second_line[56:71] 
				### now read the first line of that 
				sector_prefixes[sector], sector_suffixes[sector] = sector_prefix, sector_suffix
				#print('sector_prefix = ', sector_prefix)
				#print('len(sector_prefix) = ', len(sector_prefix))




		#### NEW MARCH 30th 2022 -- check for a record of known unobserved sectors -- so you only have to do this once!
		unobserved_sectors_filename = 'unobserved_sectors.txt'
		if os.path.exists(download_directory+'/'+unobserved_sectors_filename):
			unobserved_sectors_file = open(download_directory+'/'+unobserved_sectors_filename, mode='r')
			unobserved_sectors = []
			for nline,line in enumerate(unobserved_sectors_file):
				unobserved_sectors.append(int(line.split()[1]))
			unobserved_sectors_file.close()
		else:
			unobserved_sectors = []



		for sector in sector_numbers:
			if sector in empty_curlscript_sectors:
				continue


			if os.path.exists(central_data_dir+'/TESS_lightcurves') == False:
				os.system('mkdir '+central_data_dir+'/TESS_lightcurves')
			#download_directory = central_data_dir+'/TESS_lightcurves/TIC'+str(ticnum)

			lcdownload_name = 'TIC'+ticnum+'_sector'+str(sector)+'-s_lc.fits'
			if os.path.exists(download_directory+'/'+lcdownload_name):
				print('file already exists for sector '+str(sector))	
			else:
				if sector in unobserved_sectors:
					print('target is known not to have been observed in sector '+str(sector))
					continue

				else:
					print('attempting to download: ', lcdownload_name)
					os.system('curl  -s -C - -L -o '+download_directory+'/'+lcdownload_name+' https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/'+sector_prefixes[sector]+query_num+sector_suffixes[sector])

			try:
				lcfile = pyfits.open(download_directory+'/'+lcdownload_name)
			
			except:
				print(lcdownload_name+' was an empty data file. removing.')
				os.system('rm -rf '+download_directory+'/'+lcdownload_name)
				#### add this to unobserved_sectors_file
				unobserved_sectors_file = open(download_directory+'/'+unobserved_sectors_filename, mode='a')
				unobserved_sectors_file.write('sector '+str(sector)+'\n')
				print("sector added to 'unobserved_sectors.txt' file.")
				unobserved_sectors_file.close()

				print(' ')
				continue

			lcfiles.append(lcfile)
			lcdata = lcfile[1].data
			lctimes = np.array(lcdata['TIME'])
			if lc_format == 'pdc':
				lcfluxes = np.array(lcdata['PDCSAP_FLUX'])
				lcerrors = np.array(lcdata['PDCSAP_FLUX_ERR'])
			elif lc_format == 'sap':
				lcfluxes = np.array(lcdata['SAP_FLUX'])
				lcerrors = np.array(lcdata['SAP_FLUX_ERR'])
			lcflags = np.array(lcdata['QUALITY'])
			sector = lcfile[0].header['SECTOR']

			all_times.append(lctimes)
			all_fluxes.append(lcfluxes)
			all_errors.append(lcerrors)
			all_flags.append(lcflags)
			sectors.append(sector)

			if delete_fits == 'y':
				#os.system('rm -rf '+moonpydir+'/TESS_lcs/'+lcdownload_name)
				os.system('rm -rf '+download_directory+'/'+lcdownload_name)


	except:
		traceback.print_exc()
		print('EXCEPTION RAISED FOR tess_target_download().')
		time.sleep(2)

		obsTable = Observations.query_object(targID, radius='0.001 deg')
		TESS_idxs = np.where(np.array(obsTable['obs_collection']) == 'TESS')[0]
		minTESSidx, maxTESSidx = np.nanmin(TESS_idxs), np.nanmax(TESS_idxs)+1
		dataproducts = Observations.get_product_list(obsTable[minTESSidx:maxTESSidx])
		timeseries_idxs = np.where(np.array(dataproducts['dataproduct_type']) == 'timeseries')[0]
		obsids = np.array(dataproducts)['obsID'][timeseries_idxs]

		for obsid in np.unique(obsids):
			print ("obsid = ", obsid)
			dataproductsbyID = Observations.get_product_list(obsid)
			manifest = Observations.download_products(dataproductsbyID, download_dir=moonpydir+'/TESS_lcs', dataproduct_type='timeseries', extension='lc.fits', mrp_only=True)
			

			for nmanfile,manfile in enumerate(manifest):
				manfilepath = manfile[0]
				if "_lc.fits" in manfilepath:
					print('found the light curve!')
					### this is the only one you want to save!
					lcpath = manfilepath
					print("lcpath = ", lcpath)

					### open the file, grab the data!
					lcfile = pyfits.open(lcpath)
					lcfiles.append(lcfile)
					lcdata = lcfile[1].data
					lctimes = np.array(lcdata['TIME'])
					if lc_format == 'pdc':
						lcfluxes = np.array(lcdata['PDCSAP_FLUX'])
						lcerrors = np.array(lcdata['PDCSAP_FLUX_ERR'])
					elif lc_format == 'sap':
						lcfluxes = np.array(lcdata['SAP_FLUX'])
						lcerrors = np.array(lcdata['SAP_FLUX_ERR'])
					lcflags = np.array(lcdata['QUALITY'])
					sector = lcfile[0].header['SECTOR']

					all_times.append(lctimes)
					all_fluxes.append(lcfluxes)
					all_errors.append(lcerrors)
					all_flags.append(lcflags)
					sectors.append(sector)

					if delete_fits == 'y':
						os.system('rm '+lcpath)
					break

				else:
					pass
					#os.system('rm -rf '+manfilepath)


			print(" ")
			print(" ")

	all_times, all_fluxes, all_errors, all_flags, sectors = np.array(all_times, dtype=object), np.array(all_fluxes, dtype=object), np.array(all_errors, dtype=object), np.array(all_flags, dtype=object), np.array(sectors, dtype=object)

	return all_times, all_fluxes, all_errors, all_flags, sectors 




def eleanor_target_download(targID, sectors='all', short_cadence=False, lc_format='pdc'):
	print('calling mp_lcfind.py/eleanor_target_download().')
	import eleanor
	if sectors=='all':
		sector_array = np.array([1,2])

	tic_times, tic_sap_flux, tic_pdc_flux, tic_errors = [], [], [], []
	for sector in sector_array:
		try:
			ticstar = eleanor.Source(tic=targID, sector=sector)
			ticdata = eleanor.TargetData(ticstar, height=15, width=15, bkg_size=31, do_psf=True, do_pca=True)
			qflag0 = ticdata.quality == 0 
			tic_time, tic_raw_flux, tic_corr_flux, tic_error = ticdata.time[qlfag0], ticdata.raw_flux[qflag0], ticdata.corr_flux[qflag0], ticdata.flux_err[qflag0]

			tic_times.append(tic_time)
			tic_sap_flux.append(tic_raw_flux)
			tic_pdc_flux.append(tic_pdc_flux)
			tic_errors.append(tic_error)

		except:
			pass

	if lc_format=='pdc':
		return tic_times, tic_pdc_flux, tic_errors
	elif lc_format=='sap':
		return tic_times, tic_sap_flux, tic_errors 





def eleanor_coord_download(ra,dec, sectors='all', short_cadence=False):
	print('calling mp_lcfind.py/eleanor_coord_download().')
	print("nothing doing right now.")




def TESS_QLP_load(tic, sectors='all', clobber='n'):
	#### first download 
	TESS_QLP_download(tic, sectors=sectors, clobber=clobber)
	ticnum = tic 
	if ticnum.lower().startswith('tic'):
		ticnum = ticnum[3:]
	if ticnum.lower().startswith('-') or ticnum.lower().startswith(' '):
		ticnum = ticnum[1:]


	### now load them as you do with the other TESS light curves.
	all_times = []
	all_fluxes = []
	all_errors = []
	all_flags = []
	sectors = []
	lcfiles = []	

	if os.path.exists(central_data_dir+'/TESS_lightcurves/TIC_FFI_LCs') == False:
		os.system('mkdir '+central_data_dir+'/TESS_lightcurves/TIC_FFI_LCs')
	download_directory = central_data_dir+'/TESS_lightcurves/TIC_FFI_LCs/TIC'+str(ticnum)	
	QLP_files = os.listdir(download_directory)

	#for sector in np.arange(1,nsectors+1,1):
	for QLP_file in QLP_files:
		if '.fits' not in QLP_file:
			continue

		lcdownload_name = QLP_file 
		### find the sector 
		sector_number = lcdownload_name[19:23]
		if sector_number.startswith('0'):
			sector_number = sector_number[1:]
		sector = int(sector_number) 

		try:
			lcfile = pyfits.open(download_directory+'/'+lcdownload_name)
		except:
			os.system('rm -rf '+download_directory+'/'+lcdownload_name)
			continue

		lcfiles.append(lcfile)
		lcdata = lcfile[1].data
		lctimes = np.array(lcdata['TIME'])
		#if lc_format == 'pdc':
		lcfluxes = np.array(lcdata['SAP_FLUX'])
		lcerrors = np.array(lcdata['KSPSAP_FLUX_ERR'])
		#elif lc_format == 'sap':
		#	lcfluxes = np.array(lcdata['SAP_FLUX'])
		#	lcerrors = np.array(lcdata['SAP_FLUX_ERR'])
		lcflags = np.array(lcdata['QUALITY'])
		#sector = lcfile[0].header['SECTOR']

		all_times.append(lctimes)
		all_fluxes.append(lcfluxes)
		all_errors.append(lcerrors)
		all_flags.append(lcflags)
		sectors.append(sector)


	all_times, all_fluxes, all_errors, all_flags, sectors = np.array(all_times), np.array(all_fluxes), np.array(all_errors), np.array(all_flags), np.array(sectors)

	return all_times, all_fluxes, all_errors, all_flags, sectors 




















def TESS_QLP_download(tic, sectors='all', clobber='n'):
	if os.path.exists(central_data_dir+'/TESS_lightcurves') == False:
		os.system('mkdir '+central_data_dir+'/TESS_lightcurves')
	TESSdir = central_data_dir+'/TESS_lightcurves'
	TIC_FFI_LCs = TESSdir+'/TIC_FFI_LCs'

	if tic.lower().startswith('tic'):
		tic = tic[3:]
	if tic.lower().startswith('-') or tic.lower().startswith(' '):
		tic = tic[1:]

	TIC_filedir = TIC_FFI_LCs+'/TIC'+str(tic)

	if os.path.exists(TIC_filedir):
		nfiles_present = len(os.listdir(TIC_filedir))
		if  nfiles_present > 0:
			#### files exist already
			if clobber == 'n':
				print(str(nfiles_present)+" already present. Set clobber='y' if you think there should be more.")
				proceed_to_download = 'n'
			else:
				proceed_to_download = 'y'

		else:
			proceed_to_download = 'y'
	else:
		proceed_to_download = 'y'

	if proceed_to_download == 'y':
		if sectors == 'all':
			sector_nums = np.array([1, 14, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26])
		else:
			sector_nums = sectors 

		files_exist = '?'

		skip_sectors = np.array([])

		for sn in sector_nums:
			if sn in skip_sectors:
				continue

			print('sn = ', sn)

			sector_string = str(sn) #### will be of the form S0001 for sector 1 or S0026 for sector 26, etc.
			while len(sector_string) < 4:
				sector_string = '0'+sector_string
			assert len(sector_string) == 4
			sector_string = 's'+sector_string

			tid_string = str(tic)
			### from Huang et al, tid_string must be 16 digits, zero-padded.	
			tid_string = (16-len(tid_string))*'0'+tid_string
			assert len(tid_string) == 16
			
			txt_lcfilename = 'hlsp_qlp_tess_ffi_'+sector_string+'-'+tid_string+'_tess_v01_llc.txt'
			fits_lcfilename = 'hlsp_qlp_tess_ffi_'+sector_string+'-'+tid_string+'_tess_v01_llc.fits'

			#### here's the formatting:
			#https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/qlp/s0026/0000/0001/6554/9165/hlsp_qlp_tess_ffi_s0026-0000000165549165_tess_v01_llc.fits

			uri = 'mast:HLSP/qlp/'+sector_string+'/0000/0001/'+str(tic)[-8:-4]+'/'+str(tic)[-4:]+'/'+fits_lcfilename 
			download_URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri='+uri

			#if download_all == 'n':
			request = requests.get(download_URL)
			if request.status_code == 200:
				#### file exists!
				print('uri exists: '+str(uri))
				uri_exists = 'y'
				files_exist = 'y'	

				if (sn > 0) and (sn < 14):
					#### implies the hemisphere
					skip_sectors = np.arange(14,27,1)
					print("skipping sectors: ", skip_sectors)
				elif (sn >= 14) and (sn < 27):
					#### implies the hemisphere
					skip_sectors = np.arange(1,14,1)
					print('skipping sectors: ', skip_sectors)

			else:
				uri_exists = 'n'

			if files_exist == 'y':
				#### make sure there's a place to put the file first!
				if os.path.exists(TIC_FFI_LCs):
					pass
				else:
					os.system('mkdir '+TIC_FFI_LCs)

				if os.path.exists(TIC_FFI_LCs+'/TIC'+str(tic)):
					pass
				else:
					os.system('mkdir '+TIC_FFI_LCs+'/TIC'+str(tic))

				TIC_filepath = TIC_FFI_LCs+'/TIC'+str(tic)+'/'+fits_lcfilename
				
				try:
					#### get the file!
					urlretrieve(download_URL, TIC_filepath)	
				
				except:				
					print('URL not found for sector '+str(sn))
					continue



