from __future__ import division
import numpy as np
from moonpy import *
from astropy.io import ascii
import os
from scipy.ndimage import median_filter
from scipy.interpolate import interp1d
import time
from mp_detrend import *
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 

### FROM ASTROPY
eq_RSun = R_sun.value ### meters
eq_RJup = R_jup.value ### meters
eq_Rearth = R_earth.value ### meters
MEarth = M_earth.value ### kg
MJup = M_jup.value ### kg
MSun = M_sun.value ### kg 
AU_meters = 1.496e11 


#scrub_cnnlc_dir = input("DO YOU WANT TO SCRUB THE CNN LC DIRECTORY Y/N: ")
#if scrub_cnnlc_dir == 'y':
#	print("this will skip generating new cnn lc files and ONLY SCRUB WHAT'S ALREADY BEEN GENERATED.")
#	scrub_cnnlc_dir = input('ARE YOU SURE YOU WANT TO SCRUB THE CNN LC DIRECTORY? Y/N: ')

show_plots = input('Show plots? y/n: ')
send_to_umbriel = input('Send light curves to umbriel? y/n: ')

generate_mask_test_lcs = input('Do you want to generate the mask test light curves? y/n: ')
if generate_mask_test_lcs == 'n':
	use_randomized_transits = input('Do you want to use the randomized transit injections? y/n: ')
else:
	use_randomized_tranits = 'n'
#replace_contaminants = input('Replace contaminants with a median filter interpolation? y/n: ')
#if replace_contaminants == 'y':
#	run_all_lcs = input("Do you want to detrend 'a'll light curves, or just the 'c'ontaminated ones? ")

moonpydir = '/Users/hal9000/Documents/Software/MoonPy'
#if replace_contaminants == 'y':
#	cnnlcdir = moonpydir+'/cnn_lcs_contams_removed'
#elif replace_contaminants == 'n':
#	cnnlcdir = moonpydir+'/cnn_lcs'
if use_randomized_transits == 'n':
	if generate_mask_test_lcs == 'y':
		cnnlcdir = moonpydir+'/fpmask_test_cnn_lcs'
	else:
		cnnlcdir = moonpydir+'/fptest_cnn_lcs' ### destination for the CNN segments
	lcdir = '/Users/hal9000/Documents/Software/MoonPy/fptest_planet_injections' ### source of your light curves
elif use_randomized_transits == 'y':
	cnnlcdir = moonpydir+'/fptest_cnn_lcs_random_loc' ### destination for the CNN segments
	lcdir = '/Users/hal9000/Documents/Software/MoonPy/fptest_planet_injections_random_loc' ### source of your light curves
all_lcfiles = os.listdir(lcdir) ### every light curve from which you're generating CNN segments.

### this script will generate CNN-ready light curves for every single KOI!
cumkois = ascii.read(moonpydir+'/cumkois.txt')

all_kois = np.array(cumkois['kepoi_name'])
kepois = np.array(cumkois['kepoi_name'])
disposition = np.array(cumkois['koi_disposition'])
periods = np.array(cumkois['koi_period']).astype(float) ### days
smas = np.array(cumkois['koi_sma']).astype(float) ### AU 
smas_meters = smas * AU_meters
tau0s = np.array(cumkois['koi_time0bk']).astype(float) ### BKJD 
impacts = np.array(cumkois['koi_impact']).astype(float) ### dimensionless
eccens = np.array(cumkois['koi_eccen']).astype(float) ### dimensionless
longps = np.array(cumkois['koi_longp']).astype(float)
incs = np.array(cumkois['koi_incl']).astype(float)
durations_hours = np.array(cumkois['koi_duration']).astype(float) ### hours 
durations_days = durations_hours / 24
rprstars = np.array(cumkois['koi_ror']).astype(float) ### dimensionless
rp_rearths = np.array(cumkois['koi_prad']).astype(float)
rp_meters = rp_rearths * eq_Rearth
rstar_meters = rp_meters / rprstars ### R*[m] = (R*/Rp) * Rp[m] = Rp[m] / (Rp/R*)
lda1s = np.array(cumkois['koi_ldm_coeff1']).astype(float) 
lda2s = np.array(cumkois['koi_ldm_coeff2']).astype(float)
a_over_rstars = smas_meters / rstar_meters

already_generated = np.array(os.listdir(cnnlcdir))
already_examined = []

try:
	### this will run if these files... if they don't, they'll be created below!
	if use_randomized_transits == 'n':
		if generate_mask_test_lcs == 'y':
			already_examined_file = open(moonpydir+'/fpmask_test_already_examined.txt', mode='r')
		else:
			already_examined_file = open(moonpydir+'/fptest_already_examined.txt', mode='r')
	elif use_randomized_transits == 'y':
		already_examined_file = open(moonpydir+'/fptest_random_loc_already_examined.txt', mode='r')
	
	for naefl, aefl in enumerate(already_examined_file):
		already_examined.append(aefl)
	already_examined_file.close()
except:
	print("already_examined file has not yet been created.")
	pass 


print('you have already generated cnnlc segment files for ', len(already_generated), ' KOIs.')
print('there are a total of ', len(all_kois)-len(already_generated), 'KOIs left to process.')
print('there are a total of ', len(all_lcfiles), 'light curves left to process.')
#print('you have ', len(all_lcfiles) - len(already_generated), "left to process.")
continue_query = input('Do you wish to continue? y/n: ')
if continue_query != 'y':
	raise Exception('you opted not to continue.')

#nsegs_generated = 0
for nkoi, koi in enumerate(all_kois):

	koi_period = periods[nkoi]
	if koi_period < 5: 
		skip = 'y'
		print('this KOI has a period less than 5 days. SKIPPING.')
		continue 

	try:
		skip_koi = 'n'
		lcname = koi
		while lcname.startswith('K') or lcname.startswith('0'):
			lcname = lcname[1:]
		lcname = "KOI-"+lcname

		koi_cnnlcdir = cnnlcdir+'/'+lcname #### the directory where you'll put this CNN when it's generated!

		print("processing "+lcname)
		#bad_kois = ['3540', '3566', '3600']
		#for bk in bad_kois:
		#	if bk in lcname:
		#		skip_koi = 'y'
		#		break
		#if skip_koi == 'y':
		#	continue

		### check whether this KOI has already had CNN files generated.
		skip = 'n'
		if (lcname in already_examined) or ((lcname+'\n') in already_examined):
			skip = 'y'
			print('already examined.')
			print(' ')
			continue


		else:
			### NOT IN the ALREADY EXAMINED FILE -- TEST IF IT'S IN ALREADY GENERATED LIST.
			for ag in already_generated:
				if ag.startswith(lcname):
					if (lcname not in already_examined) and ((lcname+'\n') not in already_examined):
						### add it if it's not already there.
						#if replace_contaminants == 'y':
						#	already_examined_file = open(moonpydir+'/already_examined_contams_removed.txt', mode='a')
						#elif replace_contaminants == 'n':
						#	already_examined_file = open(moonpydir+'/already_examined.txt', mode='a')
						if use_randomized_transits == 'n':
							if generate_mask_test_lcs == 'y':
								already_examined_file = open(moonpydir+'/fpmask_test_already_examined.txt', mode='a')
							else:
								already_examined_file = open(moonpydir+'/fptest_already_examined.txt', mode='a')
						elif use_randomized_transits == 'y':
							already_examined_file = open(moonpydir+'/fptest_random_loc_already_examined.txt', mode='a')
						already_examined_file.write(lcname+'\n')
						already_examined_file.close()
					skip = 'y'
					print('already processed. Added to already_examined.txt.')
					print(' ')
					break

		if skip == 'y':
			continue

		### LOAD THE LIGHT CURVE FILE!
		if use_randomized_transits == 'n':
			if generate_mask_test_lcs == 'n':
				lcfilename = lcdir+'/'+lcname+'_kepler_fptest_lightcurve.tsv' ### this is only if you haven't randomized locations!
				lcfilenames = np.array([lcfilename])
				all_iterations.append(0)

			elif generate_mask_test_lcs == 'y':
				lcfilenames = []
				all_iterations = []
				for mask in np.arange(0,493,1):
					lcfilename = lcname+'_kepler_fpmask_test_lightcurve_mask'+str(mask)+'.tsv'
					if lcfilename in all_lcfiles:
						lcfilenames.append(lcfilename)
						all_iterations.append(mask)

		elif use_randomized_transits == 'y':
			#### in this case, there may be more than one iteration! gotta catch em all!
			#all_lcfiles = os.listdir(lcdir)

			lcfilenames = []
			all_iterations = []
			for iteration in np.arange(0,150,1): ### maximum number of iterations == maximum number of sims.
				lcfilename = lcname+'_kepler_fptest_lightcurve_iteration'+str(iteration)+'.tsv'
				if lcfilename in all_lcfiles:
					lcfilenames.append(lcfilename)
					all_iterations.append(iteration)
				else:
					break



		for iteration, lcfn in zip(all_iterations, lcfilenames):
			lcfile = open(lcdir+'/'+lcfn, mode='r')
			print('lc filename = ', lcfn)
			print('iteration = ', iteration)

			BKJD, fluxes, errors, flags, quarter, in_transit, transiter = [], [], [], [], [], [], [] ### will be arrays of arrays, to match moonpy formatting
			previous_quarter = -1
			for nline, line in enumerate(lcfile):
				linesplit = line.split()

				if nline == 0:
					try:
						new_tau0 = float(linesplit[-1]) ### this is the tau0 for the shifted transits!!!!!!!
					except:
						new_tau0 = linesplit[-1]
						while new_tau0[-1].isalpha():
							new_tau0 = new_tau0[:-1]
						new_tau0 = float(new_tau0)
					continue

				if nline == 1:
					continue
				
				else:
					this_quarter = int(float(linesplit[3]))
					
					if this_quarter > previous_quarter:
						start_new_quarter = 'y'
						#quarter.append(this_quarter)
						try:
							BKJD.append(np.array(quarter_times))
							fluxes.append(np.array(quarter_fluxes))
							errors.append(np.array(quarter_errors))
							flags.append(np.array(quarter_flags))
							in_transit.append(np.array(quarter_in_transit))
							transiter.append(np.array(quarter_transiter))
							quarter.append(this_quarter)
						except:
							pass
						previous_quarter = this_quarter
						quarter_times, quarter_fluxes, quarter_errors, quarter_flags, quarter_in_transit, quarter_transiter = [], [], [], [], [], []
					
					else:
						start_new_quarter = 'n'

					quarter_times.append(float(linesplit[0]))
					quarter_fluxes.append(float(linesplit[1]))
					quarter_errors.append(float(linesplit[2]))
					quarter_flags.append(0)
					#quarter.append(int(linesplit[3]))
					quarter_in_transit.append(linesplit[4])
					quarter_transiter.append(linesplit[5])

			BKJD = np.array(BKJD)
			fluxes = np.array(fluxes)
			errors = np.array(errors)
			flags = np.array(flags)
			quarter = np.array(quarter)
			in_transit = np.array(in_transit)
			transiter = np.array(transiter)

			koi_period = periods[nkoi]
			koi_tau0 = tau0s[nkoi]
			koi_impact = impacts[nkoi]
			koi_duration_hours = durations_hours[nkoi]
			koi_rprstar = rprstars[nkoi]
			koi_sma_AU = smas[nkoi]
			koi_rp_rearth = rp_rearths[nkoi]

			lc_dict = {}
			lc_dict['period'] = koi_period
			#lc_dict['tau0'] = koi_tau0
			lc_dict['tau0'] = new_tau0 
			lc_dict['impact'] = koi_impact
			lc_dict['duration_hours'] = koi_duration_hours
			lc_dict['rprstar'] = koi_rprstar
			lc_dict['sma_AU'] = koi_sma_AU
			lc_dict['rp_rearth'] = koi_rp_rearth 

			### now initiate a moonpy object!
			lcobject = MoonpyLC(lc_times=BKJD, lc_fluxes=fluxes, lc_errors=errors, lc_flags=flags, lc_quarters=np.unique(quarter), usr_dict=lc_dict, save_lc='n')
			lcobject.detrend(save_lc='n')


			### now with detrended light curve in hand, you can generate the CNN files.
			### don't exclude neighbors! Just flag them... maybe these will come in handy later on.
			cnnlc_path_list = lcobject.prep_for_CNN(save_lc='y', window=6, cnn_len=493, exclude_neighbors='n', flag_neighbors='n', show_plot='n', extra_path_info='iteration'+str(iteration), cnnlc_path=koi_cnnlcdir)
			nsegs_generated = len(cnnlc_path_list)

			### RENAME ALL THESE LIGHT CURVE FILES BEFORE CONTINUING!
			for ncnnlc,cnnlc in enumerate(cnnlc_path_list):
				original_path = cnnlc
				usr_planet_name = original_path[original_path.find('USR'):original_path.find('_transit')]
				new_path = original_path.replace(usr_planet_name, lcname+'_fptest')
				cnnlc_path_list[ncnnlc] = new_path 
				os.system('mv '+original_path+' '+new_path)


			if send_to_umbriel == 'y':
				remote_dir = '/home/cal/ateachey/Documents/Kepler/machine_learning_moons/fptest_cnnlc_arrays_by_planet/'+lcname
				os.system('ssh ateachey@umbriel.astro.columbia.edu mkdir '+remote_dir)
				### finally, send a copy to umbriel!
				#if replace_contaminants == 'y':
				#	remote_dir = '/home/cal/ateachey/Documents/Kepler/machine_learning_moons/cnnlc_arrays_contamfree'
				#elif replace_contaminants == 'n':
				#	remote_dir = '/home/cal/ateachey/Documents/Kepler/machine_learning_moons/cnnlc_arrays'
				

				for cnnlcfile in cnnlc_path_list:
					os.system("sftp ateachey@umbriel.astro.columbia.edu:"+remote_dir+" <<< $'put "+cnnlcfile+"'")
					print('sent to umbriel.')
				print(' ')



	except:
		traceback.print_exc()
		print('SOMETHING WENT WRONG WITH THIS KOI.')
		print('spent a lot of time trying to work out these bugs, but some of them are just too deep.')
		print("SKIPPING.")
		time.sleep(5)




