from __future__ import division
import kplr
import eleanor
import numpy as np
import astropy
from astroquery.simbad import Simbad 
import astropy.coordinates as coord 



def kplr_target_download(targID, type='koi', quarters='all', lc_format='pdc', sc=False):
	#print("nothing happening right now.")
	client = kplr.API()

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

	if coord_format=='degrees':
		nearby_objects = Simbad.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), radius='0d0m5s')
	elif coord_format=='sexagesimal':
		nearby_objects = Simbad.query_region(coord.SkyCoord(str(ra)+' '+str(dec), frame='icrs'), radius='0d0m'+str(search_radius)+'s')

	best_hit = nearby_objects[0][0].decode('UTF-8')
	if best_hit.startswith('KOI'):
		object_number = float(best_hit[4:]) ### strips off "KOI-"
		kplr_target_download(object_number, type='koi', quarters=quarters, lc_format=lc_format, sc=sc)

	elif best_hit.startswith('Kepler'):
		object_number = str(int(best_hit[7:])) ### strips off "Kepler-"
		object_name = object_number+'b'
		kplr_target_download(object_name, type='planet', quarters=quarters, lc_format=lc_format, sc=sc)






def eleanor_target_download(targID, sectors='all', sc=False):
	print("nothing happening right now.")

def eleanor_coord_download(ra,dec, sectors='all', sc=False):
	print("Nothing happening right now.")

