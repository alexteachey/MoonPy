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

moonpydir = os.getcwd()

central_data_dir = '/data/tethys/Documents/Central_Data/'


#### THIS IS A MAJOR REWORKING OF THE ORIGINAL mp_lcfind.py (eventually mp_lcfind_deprecated.py).
#### THE KEY IS, The Kplr and k2plr packages are old, not available through conda, and cause all kinds of problems.
#### So this script will get around them entirely, doing it all under the hood, the way you want it.
#### (Note that you'll want to continue structuring the light curves as kplr does it.)
#### 


#### WILL NEED TO FOLLOW THESE INSTRUCTIONS: https://astroquery.readthedocs.io/en/latest/mast/mast.html



#####################################

# XXXXXXXXXX ALIAS FINDING XXXXXXXX #

#####################################


def find_KIC_alias(target_name):
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


def find_EPIC_alias(target_name):
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



def kepler_URL_generator(KIC):
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

	download_directory = central_data_dir+'Kepler_lightcurves/KIC'+str(query_format_number)
	if os.path.exists(download_directory):
		pass
	else:
		os.system('mkdir '+download_directory)

	wget_command = "wget -q -nH --cut-dirs=6 -r -l0 -c -N -np -R --reject tar --accept fits 'index*' -erobots=off "+final_URL+" -P "+download_directory+"/"

	return final_URL, wget_command, download_directory



def k2_URL_generator(EPIC):
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

	download_directory = central_data_dir+'K2_lightcurves/EPIC'+str(query_format_number)
	if os.path.exists(download_directory):
		pass
	else:
		os.system('mkdir '+download_directory)


	#### since you don't know which campaign it's in, generate a list of final_URLs, a list of wget_commands,
	final_URLs = []
	wget_commands = []


	for campaign_number in np.arange(0,19,1):
		print('campaign # ', campaign_number)
		if campaign_number < 10:
			final_URL = 'https://archive.stsci.edu/pub/k2/lightcurves/c'+str(campaign_number)+'/'+XXXX+'00000/'+YY+'000/ktwo'+str(query_format_number)+'-c0'+str(campaign_number)+'_llc.fits'
		else:
			final_URL = 'https://archive.stsci.edu/pub/k2/lightcurves/c'+str(campaign_number)+'/'+XXXX+'00000/'+YY+'000/ktwo'+str(query_format_number)+'-c'+str(campaign_number)+'_llc.fits'
		final_URLs.append(final_URL)

		#wget_command = "wget -q -nH --cut-dirs=6 -r -l0 -c -N -np -R 'index*' -erobots=off "+final_URL+" -P "+download_directory+"/"
		#wget_command = 'wget -nH --cut-dirs=6 -r -l0 -c -N -np -R "index*" -erobots=off '+final_URL+' -P '+download_directory+'/'
		#os.system('wget --tries=1 "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_period_err1,koi_period_err2,koi_sma,koi_sma_err1,koi_sma_err2,koi_insol,koi_insol_err1,koi_insol_err2,koi_time0bk,koi_time0bk_err1,koi_time0bk_err2,koi_impact,koi_impact_err1,koi_impact_err2,koi_duration,koi_duration_err1,koi_duration_err2,koi_eccen,koi_eccen_err1,koi_eccen_err2,koi_longp,koi_longp_err1,koi_longp_err2,koi_ror,koi_ror_err1,koi_ror_err2,koi_incl,koi_incl_err1,koi_incl_err2,koi_prad,koi_prad_err1,koi_prad_err2,koi_ldm_coeff2,koi_ldm_coeff1,ra,dec&order=kepoi_name&format=ascii" -O "'+kep_mast_address+'"')
		wget_command = 'wget --tries=1 -N "'+final_URL+'" -P '+download_directory+'/'
		wget_commands.append(wget_command)


	return final_URLs, wget_commands, download_directory










#####################################

# XXXXXXXXXX DOWNLOADING XXXXXXXX #

#####################################




def kepler_fits_download(target_name):
	### first extract the KIC number
	try:
		KIC_name = find_KIC_alias(target_name)

		KIC_URL, KIC_wget, KIC_download_dir = kepler_URL_generator(KIC_name)

		print('wgetting KIC light curves...')
		os.system(KIC_wget)
		print('done.')
	except:
		traceback.print_exc()
		raise Exception("something went wrong in 'kepler_fits_download()' (see traceback).")


def k2_fits_download(target_name):
	try:
		EPIC_name = find_EPIC_alias(target_name)

		EPIC_URLs, EPIC_wgets, EPIC_download_dir = k2_URL_generator(EPIC_name) ### EPIC_URLs and EPIC_wgets are LISTS!!!!

		print('wgetting EPIC light curves...')
		for EPIC_URL, EPIC_wget in zip(EPIC_URLs, EPIC_wgets):
			print('looking for ', EPIC_URL)
			os.system(EPIC_wget)
		print('done.')
	except:
		traceback.print_exc()
		raise Exception("something went wrong in 'k2_fits_download()' (see traceback).")



#####################################

# XXXXXXXXXX UNPACK FITS XXXXXXXX #

#####################################



def kepler_unpack_fits(target_name, sc=False):
	#### first thing you have to do is make sure the damn thing exists!
	KIC_directory = kepler_URL_generator(find_KIC_alias(target_name))[2] 
	try:
		#### use kepler_fits_download() to grab the KIC_download_dir and unpack these fits.
		#KIC_directory = kepler_fits_download(target_name, download=download)
		KIC_directory_files = os.listdir(KIC_directory)

		### start a dictionary with quarters, times, fluxes, errors, flags
		kic_quarters_dict = {}

		KIC_fits_files = []

		### if sc=True, you need to check whether there is a slc version!
		### if sc=False, reject slcs!
		reject_wrong_cadence_files = []

		for i,ki in enumerate(KIC_directory_files):
			for j,kj in enumerate(KIC_directory_files):
				if (i != j) and (ki[:-8] == kj[:-8]) and (ki[-8:] == 'slc.fits') and (kj[-8:] == 'llc.fits'):
					if sc == True:
						reject_wrong_cadence_files.append(kj)
					else:
						reject_wrong_cadence_files.append(ki)


		for Kdf in KIC_directory_files:

			if ('.fits' in Kdf) and (Kdf not in reject_wrong_cadence_files):
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





def kplr_target_download(targID, targtype='koi', quarters='all', lc_format='pdc', telescope='kepler', sc=False):
	#### using the functions developed above
	### first, try unpacking without downloading.
		
	### FOR MOONPY COMPATIBILITY, THIS REQUIRES SUPPORT FOR BOTH KEPLER AND K2 LIGHT CURVES.

	if telescope.lower() == 'kepler':
		kepler_fits_download(targID)
		try:
			kepler_lc_dictionary = kepler_unpack_fits(targID, sc=sc)
		except:
			print('first except triggered.')
			time.sleep(10)
			try:
				print('lc may not have been downloaded. Attempting to download...')
				kepler_fits_download(targID)
				### after it's been downloaded, you can try to unpack again.
				kepler_lc_dictionary = kepler_unpack_fits(targID, sc=sc)
			except:
				print('second except triggered.')
				time.sleep(10)
				traceback.print_exc()
				raise Exception("Something went wrong in 'kplr_target_download()' (see traceback).")

	elif telescope.lower() == 'k2':
		k2_fits_download(targID)
		try:
			kepler_lc_dictionary = k2_unpack_fits(targID)
		except:
			print('first except triggered.')
			time.sleep(10)
			try:
				print('lc may not have been downloaded. Attempting to download...')
				k2_fits_download(targID)
				### after it's been downloaded, you can try to unpack again.
				kepler_lc_dictionary = k2_unpack_fits(targID)
			except:
				print('second except triggered.')
				time.sleep(10)
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





def kplr_coord_download(ra, dec, coord_format='degrees', quarters='all', search_radius=5, lc_format='pdc', sc=False):
	### find the object in Simbad using it's coordinates, and call kplr_target_download

	### try to interpret the coordinate_format 
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

	if (best_hit.startswith('KOI') == False) and (best_hit.startswith('Kepler') == False) and (best_hit.startswith("KIC") == False):
		### look for aliases
		alias_query = Simbad.query_objectids(str(best_hit))['ID']
		search_idx = 0
		while search_idx < 10: ### maximum tries
			try:
				alias = alias_query[search_idx]
				if alias.startswith('KOI') or alias.statswith('Kepler') or alias.startswith("KIC"):
					best_hit = alias
					break
				else:
					search_idx += 1
			except:
				break 

		print("alias = ", alias)
		best_hit = alias

	if best_hit.startswith('KOI'):
		object_number = float(best_hit[4:]) ### strips off "KOI-"
		### it's valuable to keep the quarters separated like this, because they should be detrended separately!
		targtype = 'koi'
		object_name = "KOI-"+str(object_number)

	elif best_hit.startswith('Kepler'):
		object_number = str(int(best_hit[7:]))+'b' ### strips off "Kepler-", adds b because we don't care which planet it is.
		#object_number = object_number+'b'
		targtype = 'planet'
		object_name = 'Kepler-'+str(object_number)

	elif best_hit.startswith("KIC"):
		object_number = int(best_hit[4:])
		targtype = 'kic'
		object_name = 'KIC'+str(object_number)

	print("object_number = ", object_number)


	if lc_format == 'pdc':
		kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, sc=sc)
	elif lc_format == 'sap':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, sc=sc)
	elif lc_format == 'both':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, sc=sc)


	### it's valuable to keep the quarters separated like this, because they should be detrended separately!
	if lc_format == 'pdc':
		return kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'sap':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'both':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name







def tess_coord_download(ra, dec, coord_format='degrees', quarters='all', search_radius=5, lc_format='pdc', sc=False):
	### find the object in Simbad using it's coordinates, and call kplr_target_download

	print('ra,dec = ', ra, dec)

	### try to interpret the coordinate_format 
	if 'h' in str(ra) or ':' in str(ra):
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

	if (best_hit.startswith('TOI') == False) and (best_hit.startswith('toi') == False) and (best_hit.startswith('TIC') == False) and (best_hit.startswith('tic') == False):
		### look for aliases if TOI or TIC number was not listed to start.
		alias_query = Simbad.query_objectids(str(best_hit))['ID']
		search_idx = 0
		while search_idx < 10: ### maximum tries
			try:
				alias = alias_query[search_idx]
				print('potential alias: ', alias)
				if alias.startswith('TOI') or alias.statswith('TIC'):
					best_hit = alias
					break
				else:
					search_idx += 1
			except:
				break 

		print("alias = ", alias)
		best_hit = alias

	if best_hit.startswith('TOI'):
		object_number = float(best_hit[4:]) ### strips off "KOI-"
		### it's valuable to keep the quarters separated like this, because they should be detrended separately!
		targtype = 'toi'
		object_name = "TOI-"+str(object_number)


	elif best_hit.startswith("TIC"):
		object_number = int(best_hit[4:])
		targtype = 'tic'
		object_name = 'TIC'+str(object_number)

	try:
		print("object_number = ", object_number)
	except:
		pass


	if lc_format == 'pdc':
		kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = tess_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, sc=sc)
	elif lc_format == 'sap':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters = tess_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, sc=sc)
	elif lc_format == 'both':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = tess_target_download(object_number, targtype=targtype, quarters=quarters, lc_format=lc_format, sc=sc)


	### it's valuable to keep the quarters separated like this, because they should be detrended separately!
	if lc_format == 'pdc':
		return kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'sap':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'both':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name



def tess_target_download(targID, sectors='all', sc=True, lc_format='pdc', delete_fits='n'):
	### this function interfaces with MASS to download light curves based on the TIC #.
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
		if targID.startswith('TIC'):
			ticnum = str(targID)[3:]
			if ticnum.startswith(' '):
				ticnum = str(ticnum)[1:]
		else:
			ticnum = str(targID) ### you've already handled the TOI already!
			if ticnum.startswith(' '):
				ticnum = ticnum[1:]

		### prepare for the query
		query_num = ticnum
		while len(query_num) < 16:
			query_num = '0'+query_num
		assert len(query_num) == 16

		sector_prefixes, sector_suffixes = {}, {}

		### these can be found at archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html,
		### in the tesscurl_sector_NN_lc.sh files.

		nsectors = 30
		for sector in np.arange(1,nsectors,1):
			### get the curl script... then extract the prefixes and suffixes from the first line.
			try:
				sector_curl_URL = 'http://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_'+str(sector)+'_lc.sh'
				os.system('wget --tries=1 -N "'+sector_curl_URL+'" -O '+moonpydir+'/sector'+str(sector)+"_curlscript.txt")
				curltxt = open(moonpydir+'/sector'+str(sector)+'_curlscript.txt', mode='r')
				first_line = curltxt.readline()
				second_line = curltxt.readline()
				sector_prefix = second_line[16:40]
				sector_suffix = second_line[56:71] 
				### now read the first line of that 
				sector_prefixes[sector], sector_suffixes[sector] = sector_prefix, sector_suffix
				print("sector_prefix, sector_suffix = ", sector_prefix, sector_suffix)
				#nsectors += 1
			except:
				traceback.print_exc()
				break
		print('nsectors = ', nsectors)


		"""
		sector_prefixes[1], sector_suffixes[1] = 'tess2018206045859-s0001-', '-0120-s_lc.fits'
		sector_prefixes[2], sector_suffixes[2] = 'tess2018234235059-s002-', '-0121-s_lc.fits'
		sector_prefixes[3], sector_suffixes[3] = 'tess2018263035959-s0003-', '-0123-s_lc.fits'
		sector_prefixes[4], sector_suffixes[4] = 'tess2018292075959-s0004-', '-0124-s_lc.fits'
		sector_prefixes[5], sector_suffixes[5] = 'tess2018319095959-s0005-', '-0125-s_lc.fits'
		sector_prefixes[6], sector_suffixes[6] = 'tess2018349182459-s0006-', '-0126-s_lc.fits'
		sector_prefixes[7], sector_suffixes[7] = 'tess2019006130736-s0007-', '-0131-s_lc.fits'
		sector_prefixes[8], sector_suffixes[8] = 'tess2019032160000-s0008-', '-0136-s_lc.fits'
		sector_prefixes[9], sector_suffixes[9] = 'tess2019058134432-s0009-', '-0139-s_lc.fits'
		sector_prefixes[10], sector_suffixes[10] = 'tess2019085135100-s0010-', '-0140-s_lc.fits'
		sector_prefixes[11], sector_suffixes[11] = 'tess2019112060037-s0011-', '-0143-s_lc.fits'
		sector_prefixes[12], sector_suffixes[12] = 'tess2019140104343-s0012-', '-0144-s_lc.fits'
		sector_prefixes[13], sector_suffixes[13] = 'tess2019169103026-s0013-', '-0146-s_lc.fits'
		sector_prefixes[14], sector_suffixes[14] = 'tess2019198215352-s0014-', '-0150-s_lc.fits'
		sector_prefixes[15], sector_suffixes[15] = 'tess2019226182529-s0015-', '-0151-s_lc.fits'
		sector
		nsectors = 28
		"""

		for sector in np.arange(1,nsectors+1,1):
			lcdownload_name = 'TIC'+ticnum+'_sector'+str(sector)+'-s_lc.fits'
			os.system('curl -C - -L -o '+moonpydir+'/TESS_lcs/'+lcdownload_name+' https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/'+sector_prefixes[sector]+query_num+sector_suffixes[sector])
			print('downloading the light curve for '+str(targID)+' in sector ', sector)

			try:
				lcfile = pyfits.open(moonpydir+'/TESS_lcs/'+lcdownload_name)
			except:
				os.system('rm -rf '+moonpydir+'/TESS_lcs/'+lcdownload_name)
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
				os.system('rm -rf '+moonpydir+'/TESS_lcs/'+lcdownload_name)


	except:
		traceback.print_exc()
		time.sleep(60)

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

	all_times, all_fluxes, all_errors, all_flags, sectors = np.array(all_times), np.array(all_fluxes), np.array(all_errors), np.array(all_flags), np.array(sectors)

	return all_times, all_fluxes, all_errors, all_flags, sectors 




def eleanor_target_download(targID, sectors='all', sc=False, lc_format='pdc'):
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





def eleanor_coord_download(ra,dec, sectors='all', sc=False):
	print("nothing doing right now.")

