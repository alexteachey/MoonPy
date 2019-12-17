from __future__ import division
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt 
import numpy as np
import traceback
import os 
from matplotlib import animation 
from mp_detrend import *
from astroquery.simbad import Simbad 


destination_dir = '/Users/hal9000/Documents/Software/MoonPy/TPFs'
#k1654_tpf = 'kplr008410697-2010355172524_lpd-targ.fits'
#tpfdir = '/Users/hal9000/Documents/Projects/k1654/008410697/'


def find_kic(targetname):
	target_aliases = []
	alias_search_results = Simbad.query_objectids(targetname)
	for alidx in np.arange(0,np.array(alias_search_results).shape[0],1):
		target_alias = alias_search_results[alidx][0]
		target_aliases.append(target_alias)

	for ta in target_aliases:
		if ta.startswith('KIC'):
			return ta




def quadsum(vals):
	return np.sqrt(np.nansum(vals**2))


def lookup_epochs(quarter, cadence):
	LONG_QUARTER_PREFIXES = {'0':['2009131105131'], '1':['2009166043257'], '2':['2009259160929'], '3':['2009350155506'],'4':['2010078095331', '2010009091648'],
								'5':['2010174085026'],'6':['2010265121752'], '7':['2010355172524'],'8':['2011073133259'], '9':['2011177032512'], '10':['2011271113734'],
								'11':['2012004120508'], '12':['2012088054726'], '13':['2012179063303'], '14':['2012277125453'], '15':['2013011073258'], '16':['2013098041711'], '17':['2013131215648']}
    
	SHORT_QUARTER_PREFIXES = {'0':['2009131110544'], '1':['2009166044711'], '2':['2009201121230', '2009231120729', '2009259162342'], '3':['2009291181958', '2009322144938','2009350160919'],
								'4':['2010009094841', '2010019161129', '2010049094358', '2010078100744'], '5':['2010111051353', '2010140023957', '2010174090439'],'6':['2010203174610', '2010234115140','2010265121752'],
								'7':['2010296114515', '2010326094124', '2010355172524'], '8':['2011024051157', '2011053090032', '2011073133259'], '9':['2011116030358', '2011145075126','2011177032512'], '10':['2011208035123', '2011240104155', '2011271113734'],
								'11':['2011303113607', '2011334093404', '2012004120508'], '12':['2012032013838', '2012060035710', '2012088054726'], '13':['2012121044856', '2012151031540', '2012179063303'], '14':['2012211050319', '2012242122129', '2012277125453'], 
								'15':['2012310112549', '2012341132017', '2013011073258'], '16':['2013017113907', '2013065031647', '2013098041711'], '17':['2013121191144', '2013131215648']}

	if cadence == "long":
		if quarter in LONG_QUARTER_PREFIXES:
			epoch_output = LONG_QUARTER_PREFIXES[quarter]
			return epoch_output 
		else:
			raise ValueError("*** ERROR in lookup_epochs: Quarter must be between 0 and 17.")

	elif cadence == "short":
		if quarter in SHORT_QUARTER_PREFIXES:
			epoch_output = SHORT_QUARTER_PREFIXES[quarter]
			return epoch_output 
		else:
			raise ValueError("*** ERROR in lookup_epochs: Quarter must be between 0 and 17.")
	else:
		raise ValueError("*** ERROR in lookup_epochs: Cadence must be 'long' or 'short'.")



def tpf_downloader(kic, quarters, cadence='long', clobber='n'): 


	### quarters can be a single quarter or a list.

	file_dict = {}

	if (type(quarters) == str) or (type(quarters) == float):
		quarters = [str(quarters)]
	else:
		pass

	for quarterval in quarters:
		quarter = str(quarterval)
		#### get the prefixes
		quarter_codes = lookup_epochs(quarter, cadence)

		#### takes a KIC number and figures out the website where the files are hosted, etc
		if str(kic).startswith('KIC') or str(kic).startswith('kic'):
			kicnum = str(kic)[3:]
		else:
			kicnum = kic 

		if (str(kicnum).startswith(' ')) or (str(kicnum).startswith('-')):
			kicnum = kicnum[1:]

		kic_prefix = str(0)+str(0)+str(kicnum[:2]) ### first four numbers 
		full_kicnum = str(0)+str(0)+str(kicnum)
		print ('kic_prefix = ', kic_prefix)
		print('full_kicnum = ', full_kicnum)

		tp_base_url = "http://archive.stsci.edu/missions/kepler/target_pixel_files/"

		if cadence == "long":
			cadence_str = "_lpd-targ.fits.gz"
		else:
			cadence_str = "_spd-targ.fits.gz"

		kic_tpf_url = 'https://archive.stsci.edu/pub/kepler/target_pixel_files/'+kic_prefix+'/'+full_kicnum
		print('kic_tpf_url = ', kic_tpf_url)
		### at the page above, every quarter is saved as a separate tar.gz file. ugh.


		for nquarter_code, quarter_code in enumerate(quarter_codes):
			kic_tpf_filename_gz = 'kplr'+full_kicnum+'-'+str(quarter_code)+cadence_str
			kic_tpf_filename_fits = kic_tpf_filename_gz[:-3]
			print('kic_tpf_filename = ', kic_tpf_filename_gz)

			wget_kic_address = kic_tpf_url+'/'+kic_tpf_filename_gz
			print("wget address = ", wget_kic_address)

			### see whether the file already exists!
			if (os.path.exists(destination_dir+'/'+kic_tpf_filename_fits)) and (clobber == 'n'):
				print("file already exists and clobber=='n'")
				pass

			else:
				if (os.path.exists(destination_dir+'/'+kic_tpf_filename_gz)) and (clobber == 'n'):
					pass

				else:
					os.system('wget '+wget_kic_address+' -P '+destination_dir+'/')

				### unzip it
				os.system('gunzip '+destination_dir+'/'+kic_tpf_filename_gz)

				### remove the zipped version:
				os.system('rm -rf '+destination_dir+'/'+kic_tpf_filename_gz)


			fits_file_path = destination_dir+'/'+kic_tpf_filename_fits

			if len(quarter_codes) == 1:
				file_dict[quarter] = fits_file_path 
			else:
				file_dict[quarter+'_'+nquarter_code] = fits_file_path 

	return file_dict 






def tpf_examiner(kic, quarters, find_alias='n', time_lims=None, cadence='long', clobber='n', detrend='y', mask_times=None, mask_idxs=None):
	if find_alias == 'y':
		try:
			kic = find_kic(kic)
		except:
			print('could not identify a KIC alias.')

	#### kic is flexible... can be with or without 'KIC'
	### quarters may be a single value or a list.
	### time_lims, if specified, should be a tuple of (time_min, time_max) (BKJD)

	### invoke the tpf_downloader
	print(kic)
	tpf_file_dict = tpf_downloader(kic, quarters, cadence=cadence, clobber=clobber)

	for tpf_file_key in tpf_file_dict.keys():
		print('analyzing quarter ', tpf_file_key)

		tpf_path = tpf_file_dict[tpf_file_key]

		tpf = pyfits.open(tpf_path)

		times = tpf[1].data['TIME']
		qflags = tpf[1].data['QUALITY']
		raw_counts = tpf[1].data['RAW_CNTS']
		fluxes = tpf[1].data['FLUX']
		flux_errs = tpf[1].data['FLUX_ERR']
		time_correlation = tpf[1].data['TIMECORR']
		cosmic_rays = tpf[1].data['COSMIC_RAYS']

		if type(mask_times) != type(None):
			#### find the indices of these times
			mask_idxs = np.where((times >= np.nanmin(mask_times)) & (times <= np.nanmax(mask_times)))[0]

		aperture = tpf[2].data ### 0 means not collected, 1 means collected but not part of the aperture, 3 means in the aperture.

		all_active_pixel_idxs = np.where(aperture >= 1)
		aperture_pixel_idxs = np.where(aperture == 3)

		### zero in on the interesting times
		if type(time_lims) != type(None):
			window_time_idxs = np.where((times >= time_lims[0]) & (times <= time_lims[1]))[0]
		else:
			window_time_idxs = np.arange(0,len(times),1)

		window_times = []
		all_flux_lc = []
		all_flux_error_lc = []
		all_flux_CR_lc = []
		all_flux_xcentroids = []
		all_flux_ycentroids = []

		aperture_flux_lc = []
		aperture_flux_error_lc = []
		aperture_flux_CR_lc = []
		aperture_flux_xcentroids = []
		aperture_flux_ycentroids = []

		timecorrs = []
		quality_flags = []


		### try to animate here!
		"""
		nframes = len(window_time_idxs)
		fps = 10

		#fig = plt.figure(figsize = fluxes[0].shape)
		fig = plt.figure()
		all_frames = fluxes[window_time_idxs]
		a = all_frames[0]
		ims = []

		for wtidx in window_time_idxs:
			im = plt.imshow(np.nan_to_num(fluxes[wtidx]), animated=True)
			ims.append([im])

		ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True)
		plt.show()
		"""	

		"""
		im = plt.imshow(a, interpolation='none', vmin=np.nanmin(fluxes), vmax=np.nanmax(fluxes))
		im.scatter(all_active_pixel_idxs[1], all_active_pixel_idxs[0], c='k', marker='*', s=10)
		im.scatter(aperture_pixel_idxs[1], aperture_pixel_idxs[0], c='r', marker='X', s=15)

		def animate_func(i):
			im.set_array(all_frames[i], origin='lower')
			#im.scatter(all_active_pixel_idxs[1], all_active_pixel_idxs[0], c='k', marker='*', s=10)
			#im.scatter(aperture_pixel_idxs[1], aperture_pixel_idxs[0], c='r', marker='X', s=15)
			#im.set_title('BKJD = '+str(window_times[i]))
			return [im] 

		anim = animation.FuncAnimation(fig, animate_func, frames=nframes, interval=1000/fps, blit=True)
		plt.show()
		"""




		for nwtidx, wtidx in enumerate(window_time_idxs):
			#print ('wtidx = ', wtidx)
			#print('fluxes[wtidx].shape = ', fluxes[wtidx].shape)
			try:
				if qflags[wtidx] != 0:
					continue

				window_times.append(times[wtidx])

				### ALL FLUXES
				all_fluxes = fluxes[wtidx][all_active_pixel_idxs]
				all_flux_lc.append(np.nansum(all_fluxes))

				all_flux_errs = flux_errs[wtidx][all_active_pixel_idxs]
				#all_flux_error_lc.append(quadsum(all_flux_errs)) ### I don't think this should sum in quadrature
				all_flux_error_lc.append(np.nansum(all_flux_errs))

				all_flux_CRs = cosmic_rays[wtidx][all_active_pixel_idxs]
				all_flux_CR_lc.append(np.nansum(all_flux_CRs))

				### APERTURE FLUXES
				aperture_fluxes = fluxes[wtidx][aperture_pixel_idxs]
				aperture_flux_lc.append(np.nansum(aperture_fluxes))

				aperture_flux_errs = flux_errs[wtidx][aperture_pixel_idxs]
				aperture_flux_error_lc.append(np.nansum(aperture_flux_errs))

				aperture_flux_CRs = cosmic_rays[wtidx][aperture_pixel_idxs]
				aperture_flux_CR_lc.append(np.nansum(aperture_flux_CRs))



				#### COMPUTE THE X AND Y CENTROIDS!
				all_flux_centx = np.nansum(all_active_pixel_idxs[0]*fluxes[wtidx][all_active_pixel_idxs]) / np.nansum(fluxes[wtidx][all_active_pixel_idxs])
				all_flux_centy = np.nansum(all_active_pixel_idxs[1]*fluxes[wtidx][all_active_pixel_idxs]) / np.nansum(fluxes[wtidx][all_active_pixel_idxs])

				all_flux_xcentroids.append(all_flux_centx)
				all_flux_ycentroids.append(all_flux_centy)

				aperture_flux_centx = np.nansum(aperture_pixel_idxs[0]*fluxes[wtidx][aperture_pixel_idxs]) / np.nansum(fluxes[wtidx][aperture_pixel_idxs])
				aperture_flux_centy = np.nansum(aperture_pixel_idxs[1]*fluxes[wtidx][aperture_pixel_idxs]) / np.nansum(fluxes[wtidx][aperture_pixel_idxs])

				aperture_flux_xcentroids.append(aperture_flux_centx)
				aperture_flux_ycentroids.append(aperture_flux_centy)

				### make the aperture flux stack
				if nwtidx == 0:
					aperture_flux_stack = aperture_fluxes
					aperture_error_stack = aperture_flux_errs
				else:
					aperture_flux_stack = np.vstack((aperture_flux_stack, aperture_fluxes))
					aperture_error_stack = np.vstack((aperture_error_stack, aperture_flux_errs))


				### 1-D arrays
				timecorrs.append(time_correlation[wtidx])
				quality_flags.append(qflags[wtidx])


				### as a test, let's generate the plot at each time step!
				"""
				fig, ax = plt.subplots(1)
				ax.imshow(fluxes[wtidx], origin='lower')
				ax.scatter(all_active_pixel_idxs[1], all_active_pixel_idxs[0], c='k', marker='*', s=10)
				ax.scatter(aperture_pixel_idxs[1], aperture_pixel_idxs[0], c='r', marker='X', s=15)
				plt.show()
				"""

			except:
				traceback.print_exc()

		### convert to arrays
		window_times = np.array(window_times)
		all_flux_lc = np.array(all_flux_lc)
		all_flux_error_lc = np.array(all_flux_error_lc)
		all_flux_CR_lc = np.array(all_flux_CR_lc)

		aperture_flux_lc = np.array(aperture_flux_lc)
		aperture_flux_error_lc = np.array(aperture_flux_error_lc)
		aperture_flux_CR_lc = np.array(aperture_flux_CR_lc)

		all_minus_aperture_flux_lc = all_flux_lc - aperture_flux_lc 
		all_minus_aperture_flux_err_lc = all_flux_error_lc - aperture_flux_error_lc
		all_minus_aperture_CR_lc = all_flux_CR_lc - aperture_flux_CR_lc 

		timecorrs = np.array(timecorrs)
		quality_flags = np.array(quality_flags)



		### pixel fluxes
		fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
		ax1.scatter(window_times, all_flux_lc, s=10, color='LightCoral', zorder=1, label='all pixels')
		#ax1.errorbar(window_times, all_flux_lc, yerr=all_flux_error_lc, fmt='none', ecolor='k', zorder=0)
		#ax1.set_title('all recorded pixels')
		ax1.set_title('Quarter '+tpf_file_key)
		ax1.legend()
		#plt.legend()

		ax2.scatter(window_times, aperture_flux_lc, s=10, color='DodgerBlue', zorder=1, label='aperture')
		#ax2.errorbar(window_times, aperture_flux_lc, yerr=aperture_flux_error_lc, fmt='none', ecolor='k', zorder=0)
		#ax2.set_title('aperture pixels')
		ax2.legend()

		ax3.scatter(window_times, all_minus_aperture_flux_lc, color='SeaGreen', s=10, label='all - aperture')
		#ax3.set_title('all - aperture')
		ax3.set_xlabel('BKJD')
		ax3.legend()
		plt.show()


		### plot the flux stack
		for i in np.arange(0,aperture_flux_stack.shape[1], 1):
			### identify the pixel value here:
			ypix, xpix = aperture_pixel_idxs[0][i], aperture_pixel_idxs[1][i] 

			### plot the light curve for just this pixel!
			### normalize the pixel light curve.
			#normed_pixel_lc = (aperture_flux_stack.T[i] - np.nanmin(aperture_flux_stack.T[i])) / (np.nanmax(aperture_flux_stack.T[i]) - np.nanmin(aperture_flux_stack.T[i]))
			normed_pixel_lc = aperture_flux_stack.T[i] - np.nanmedian(aperture_flux_stack.T[i])

			plt.plot(window_times, normed_pixel_lc, label='x,y = '+str(xpix)+','+str(ypix))

		plt.title('Quarter '+str(tpf_file_key))
		plt.legend(loc=0)
		plt.xlabel("BKJD")
		plt.ylabel('Normalized Pixel Flux')
		plt.show()


		if detrend == 'y':
			for i in np.arange(0,aperture_flux_stack.shape[1], 1):
				ypix, xpix = aperture_pixel_idxs[0][i], aperture_pixel_idxs[1][i] 

				### detrend each individual aperture pixel individually with CoFiAM!
				pixel_flux_detrend, pixel_error_detrend = cofiam_detrend(window_times, aperture_flux_stack.T[i], aperture_error_stack.T[i], mask_idxs = mask_idxs)

				plt.plot(window_times, pixel_flux_detrend, label='x,y = '+str(xpix)+','+str(ypix))

		plt.title('Quarter '+str(tpf_file_key))
		plt.legend(loc=0)
		plt.xlabel("BKJD")
		plt.ylabel('Pixel Detrending')
		plt.show()			



		fig, (ax1, ax2) = plt.subplots(2)
		ax1.scatter(all_flux_xcentroids, all_flux_ycentroids, c=window_times, s=10)
		#ax1.colorbar('BKJD')
		#ax1.set_xlabel('x-centroid')
		ax1.set_ylabel('y-centroid')
		ax2.scatter(aperture_flux_xcentroids, aperture_flux_ycentroids, c=window_times, s=10)
		#ax2.colorbar('BKJD')
		#plt.colorbar(label="BKJD")
		ax2.set_xlabel('x-centroid')
		ax2.set_ylabel('y-centroid')
		plt.show()


		#### make the animation!


		#### plot all all pixel vs aperture pixels, color code by time
		plt.title('Quarter '+str(tpf_file_key))
		plt.scatter(all_flux_lc, aperture_flux_lc, c=window_times, s=10)
		plt.xlabel('all collected pixel fluxes')
		plt.ylabel('aperture pixel fluxes')
		plt.colorbar(label='BKJD')
		plt.show()


		#### look at the ratio vs time
		"""
		plt.title('Quarter '+str(tpf_file_key))
		plt.scatter(window_times, aperture_flux_lc / all_flux_lc, color='LightCoral', s=10)
		plt.xlabel('BKJD')
		plt.ylabel('aperture flux / all pixel flux')
		plt.show()
		"""


		#### examine cosmic rays
		fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
		ax1.set_title('Quarter '+str(tpf_file_key))
		ax1.plot(window_times, all_flux_CR_lc, color='LightCoral', label='all CRs')
		ax1.legend()
		ax2.plot(window_times, aperture_flux_CR_lc, color='DodgerBlue', label='aperture CRs')
		ax2.legend()
		ax3.plot(window_times, all_minus_aperture_CR_lc, color='SeaGreen', label='all CRs - aperture CRs')
		ax3.legend()
		ax3.set_xlabel('BKJD')
		plt.show()

		#### examine the quality flagsd

		plt.title('Quarter '+str(tpf_file_key))
		plt.plot(window_times, quality_flags, c='LightCoral')
		plt.xlabel('BKJD')
		plt.ylabel('Quality')
		plt.show()


		### examine time correlation
		plt.title('Quarter '+str(tpf_file_key))
		plt.plot(window_times, timecorrs, color='LightCoral')
		plt.xlabel('Time Correlation')
		plt.ylabel('BKJD')
		plt.show()








