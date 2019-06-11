from __future__ import division
import numpy as np
import astropy
from astroquery.simbad import Simbad 
from astroquery.mast import Observations
import astropy.coordinates as coord 
from astropy import units as u
import os
from astropy.io import fits as pyfits 

moonpydir = os.getcwd()

def kplr_target_download(targID, type='koi', quarters='all', lc_format='pdc', telescope='kepler', sc=False):
	import kplr
	import k2plr
	#print("nothing happening right now.")
	if (telescope == 'kepler') or (telescope=="Kepler"):
		client = kplr.API()
	elif (telescope == 'k2') or (telescope == "K2"):
		client = k2plr.API()

	if type == 'koi':
		kplr_obj = client.koi(targID)
	elif type == 'planet':
		kplr_obj = client.planet(str(targID))
	elif type == 'kic':
		kplr_obj = client.star(targID)

	lcs = kplr_obj.get_light_curves(short_cadence=sc)

	kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = [], [], [], [], [], [], []

	for lc in lcs:
		with lc.open() as f:
			hdu_header = f[0].header 
			hdu_data = f[1].data 

			obj_ra, obj_dec = hdu_header['RA_OBJ'], hdu_header['DEC_OBJ']
			quarter = hdu_header['QUARTER']

			kobj_times.append(hdu_data['time'])
			
			kobj_sap_fluxes.append(hdu_data['sap_flux'])
			kobj_sap_errors.append(hdu_data['sap_flux_err'])
			kobj_pdc_fluxes.append(hdu_data['pdcsap_flux'])
			kobj_pdc_errors.append(hdu_data['pdcsap_flux_err'])
			kobj_flags.append(hdu_data['sap_quality'])
			kobj_quarters.append(quarter)

	kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = np.array(kobj_times), np.array(kobj_sap_fluxes), np.array(kobj_sap_errors), np.array(kobj_pdc_fluxes), np.array(kobj_pdc_errors), np.array(kobj_flags), np.array(kobj_quarters)

	### sort them
	quarter_sort = np.argsort(kobj_quarters)

	kobj_times = kobj_times[quarter_sort]
	kobj_sap_fluxes, kobj_sap_errors = kobj_sap_fluxes[quarter_sort], kobj_sap_errors[quarter_sort]
	kobj_pdc_fluxes, kobj_pdc_errors = kobj_pdc_fluxes[quarter_sort], kobj_pdc_errors[quarter_sort]
	kobj_flags = kobj_flags[quarter_sort]
	kobj_quarters = kobj_quarters[quarter_sort]

	if quarters != 'all':
		### grab just the quarters you want
		final_quarter_idxs = []
		for nquart,quart in enumerate(kobj_quarters):
			if quart in quarters:
				final_quarter_idxs.append(nquart)
		final_quarter_idxs = np.array(final_quarter_idxs)

		kobj_times = kobj_times[final_quarter_idxs]
		kobj_sap_fluxes, kobj_sap_errors = kobj_sap_fluxes[final_quarter_idxs], kobj_sap_errors[final_quarter_idxs]
		kobj_pdc_fluxes, kobj_pdc_errors = kobj_pdc_fluxes[final_quarter_idxs], kobj_pdc_errors[final_quarter_idxs]
		kobj_flags = kobj_flags[final_quarter_idxs]
		kobj_quarters = kobj_quarters[final_quarter_idxs]


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
		kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, type=targtype, quarters=quarters, lc_format=lc_format, sc=sc)
	elif lc_format == 'sap':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, type=targtype, quarters=quarters, lc_format=lc_format, sc=sc)
	elif lc_format == 'both':
		kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters = kplr_target_download(object_number, type=targtype, quarters=quarters, lc_format=lc_format, sc=sc)


	### it's valuable to keep the quarters separated like this, because they should be detrended separately!
	if lc_format == 'pdc':
		return kobj_times, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'sap':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_flags, kobj_quarters, object_name
	elif lc_format == 'both':
		return kobj_times, kobj_sap_fluxes, kobj_sap_errors, kobj_pdc_fluxes, kobj_pdc_errors, kobj_flags, kobj_quarters, object_name




def tess_target_download(targID, sectors='all', sc=False, lc_format='pdc', delete_fits='y'):
	### this function interfaces with MASS to download light curves based on the TIC #.
	if os.path.exists(moonpydir+'/TESS_lcs'):
		pass
	else:
		os.system('mkdir '+moonpydir+'/TESS_lcs')

	obsTable = Observations.query_object(targID, radius='0.001 deg')
	TESS_idxs = np.where(np.array(obsTable['obs_collection']) == 'TESS')[0]
	minTESSidx, maxTESSidx = np.nanmin(TESS_idxs), np.nanmax(TESS_idxs)+1
	dataproducts = Observations.get_product_list(obsTable[minTESSidx:maxTESSidx])
	timeseries_idxs = np.where(np.array(dataproducts['dataproduct_type']) == 'timeseries')[0]
	obsids = np.array(dataproducts)['obsID'][timeseries_idxs]

	all_times = []
	all_fluxes = []
	all_errors = []
	all_flags = []
	sectors = []
	lcfiles = []

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

