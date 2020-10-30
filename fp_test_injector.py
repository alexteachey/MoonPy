from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas
from moonpy import *
from astropy.io import ascii
from mp_batman import *
from mp_tools import *
from scipy.signal import medfilt
from scipy.interpolate import interp1d 
import batman
from astropy.constants import G, c, M_earth, M_jup, M_sun, R_earth, R_jup, R_sun, au 
import traceback
import time
import george
from george.kernels import ExpSquaredKernel
from scipy.interpolate import interp1d 
import celerite
from celerite import terms 
from scipy.optimize import minimize
import time



### THIS CODE WILL PRODUCE A FALSE POSITIVE TEST FOR YOUR CNN VETTING.
"""
Here's the idea:
inject planets ONLY -- the same size, duration, etc as the real planet.
inject these planets into CLEAN SEGMENTS OF THE VERY SAME LIGHT CURVE -- thus, capturing the unique 
instrumental and astrophysical artifacts.
THEN we can run these through the CNN and determine how often it mis-classifies PLANET ONLY as PLANET+MOON.

Strategy:
-FOR EACH KOI, we want to 1) download the light curve (use MoonPy)
-identify real transits in the system, and mask them
-inject planets of the same size and duration, etc.
-save the light curves! (send to umbriel) -- make sure they're formatted correctly

"""



def neg_log_like(params, y, gp):
	gp.set_parameter_vector(params)
	return -gp.log_likelihood(y)


def GPfit(times, fluxes, errors, nonmask_idxs):
	#t, y, yerr = times[np.isfinite(times)], fluxes[np.isfinite(fluxes)], errors[np.isfinite(errors)]
	t, y, yerr = times[nonmask_idxs], fluxes[nonmask_idxs], errors[nonmask_idxs]
	t, y, yerr = t[np.isfinite(t)], y[np.isfinite(y)], yerr[np.isfinite(yerr)]

	Q = 1.0
	w0 = 3.0
	S0 = np.var(y) / (w0 * Q)
	bounds = dict(log_S0=(-15, 15), log_Q=(-15, 15), log_omega0=(-15, 15))
	kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0), bounds=bounds)
	#kernel += terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0) bounds=bounds)


	gp = celerite.GP(kernel, mean=np.mean(y))
	gp.compute(t, yerr)  # You always need to call compute once.
	print("Initial log likelihood: {0}".format(gp.log_likelihood(y)))

	initial_params = gp.get_parameter_vector()
	bounds = gp.get_parameter_bounds()

	r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(y, gp))
	gp.set_parameter_vector(r.x)
	print(r)

	#interpolated_times = np.linspace(np.nanmin(t), np.nanmax(t), 1000)
	#interpolated_epochs = np.linspace(np.nanmin(t), np.nanmax(t), 1000)
	#pred_mean, pred_var = gp.predict(y, interpolated_times, return_var=True)
	pred_mean, pred_var = gp.predict(y, times, return_var=True)
	pred_std = np.sqrt(pred_var)

	return pred_mean, pred_std ### should be the same dimension as times.


### important values
### constants -- THESE SHOULD BE IMPORTED FROM ASTROPY!


### FROM ASTROPY
eq_RSun = R_sun.value ### meters
eq_RJup = R_jup.value ### meters
eq_Rearth = R_earth.value ### meters
MEarth = M_earth.value ### kg
MJup = M_jup.value ### kg
MSun = M_sun.value ### kg 
AU_meters = 1.496e11 




print('This script takes real Kepler light curves, removes the original transit signals, and injects new planet signals in their place. The purpose of this is to run a false positive test for the CNN vetter.')
show_plots = input('Show plots? y/n: ')
randomize_injection_location = input('Do you want to randomize the planet injection site? y/n: ')
if randomize_injection_location == 'n':
	generate_mask_fptest_lcs = input('Do you want to generate the false positive mask light curves? ')

ntransits_per_star = int(input('How many transits do you want to generate for each light curve? '))
write_to_file = input('Do you want to write these light curves to a file? ')
overwrite_fplc = input('OVERWRITE existing files of the same name? y/n: ')
mpdir = '/Users/hal9000/Documents/Software/MoonPy'

run_individual_system = input('Do you want to run an individual system? y/n: ')

cumkois = ascii.read(mpdir+'/cumkois.txt')
kepois = np.array(cumkois['kepoi_name'])
dispositions = np.array(cumkois['koi_disposition'])
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


#run_batman(all_times, RpRstar, Rstar, bplan, Pplan, tau0, q1, q2, long_peri=0, ecc=0, Mstar=None, Mplan=None, rhostar=None, rhoplan=None, cadence_minutes=29.42, noise_ppm=None, munit='kg', runit='meters', ang_unit='radians', add_noise='n', show_plots='n', print_params='n', binned_output='n', **kwargs):

kois = []

for kepoi in kepois:
	koi = kepoi
	while koi.startswith('K') or koi.startswith('0'):
		koi = koi[1:]
	kois.append('KOI-'+str(koi))
kois = np.array(kois)


if run_individual_system == 'y':
	### HARD CODING kepler-150f in for now
	single_system = input('What is the name of the system? ')
	kois = np.array([single_system])
	if single_system.lower().startswith('koi') == False:
		iskoi = 'n'
	else:
		iskoi = 'y'



#np.random.shuffle(kois)

for nkoi, koi in enumerate(kois):

	try:
		print(koi)

		if iskoi == 'y':
			#planet_number = koi[4:koi.find('.')] ### isolates the number by itself.
			kepoi_by_itself = kepois[nkoi][:kepois[nkoi].find('.')] ### this will be something like "K000001" for KOI-1.01.
		

		### see whether there are neighbors:
		neighbor_kepois = []
		neighbor_idxs = []
		for kepoi_idx in np.arange(-5,6,1):
			try:
				if (kepoi_by_itself in kepois[nkoi+kepoi_idx]) and (kepoi_idx != 0): ### if, for example, K000001.02 exists.
					neighbor_kepoi = kepois[nkoi+kepoi_idx]
					neighbor_kepois.append(neighbor_kepoi)
					neighbor_idxs.append(nkoi+kepoi_idx)
			except:
				#print('neighbor finding failed (not necessarily bad).')
				continue
		neighbor_idxs = np.array(neighbor_idxs)

		### grab all the koi params!
		koi_period = periods[nkoi]
		koi_disposition = dispositions[nkoi]
		koi_sma = smas[nkoi]
		koi_tau0 = tau0s[nkoi]
		koi_impact = impacts[nkoi]
		koi_duration_hours = durations_hours[nkoi]
		koi_duration_days = durations_days[nkoi]
		koi_eccen = eccens[nkoi]
		koi_longp = longps[nkoi]
		koi_inc = incs[nkoi]
		koi_rprstar = rprstars[nkoi]
		koi_rp_rearth = rp_rearths[nkoi]
		koi_lda1 = lda1s[nkoi]
		koi_lda2 = lda2s[nkoi]
		koi_q1, koi_q2 = u1u2_to_q1q2(koi_lda1, koi_lda2)
		koi_rp_meters = rp_meters[nkoi]
		koi_rstar_meters = rstar_meters[nkoi]
		koi_a_over_rstar = a_over_rstars[nkoi]

		print('period = ', koi_period)
		print('impact = ', koi_impact)
		print('duration_days = ', koi_duration_days)
		print('e = ', koi_eccen)
		print('longp = ', koi_longp)
		print('inclination = ', koi_inc)

		#if float(koi_period) < 10:
		if (float(koi_period) < 5) or (float(koi_period) >= 10):
			#print('period < 10 days.')
			print('period outside the desired range. (all KOIs > 10 days already generated).')
			print(" ")
			print(" ")
			continue

		print('disposition = ', koi_disposition)
		if koi_disposition == 'FALSE POSITIVE':
			print("FALSE POSITIVE! SKIPPING.")
			#time.sleep(3)
			print(' ')
			print(' ')
			continue

		print('# of neighbors = ', len(neighbor_idxs))

		if len(neighbor_idxs) != 0:
			neighbor_periods = periods[neighbor_idxs]
			neighbor_smas = smas[neighbor_idxs]
			neighbor_tau0s = tau0s[neighbor_idxs]
			neighbor_impacts = impacts[neighbor_idxs]
			neighbor_durations_hours = durations_hours[neighbor_idxs]
			neighbor_durations_days = durations_days[neighbor_idxs]
			neighbor_eccens = eccens[neighbor_idxs]
			neighbor_longps = longps[neighbor_idxs]
			neighbor_incs = incs[neighbor_idxs]
			neighbor_rprstars = rprstars[neighbor_idxs]
			neighbor_rp_rearth = rp_rearths[neighbor_idxs]
			neighbor_lda1s = lda1s[neighbor_idxs]
			neighbor_lda2s = lda2s[neighbor_idxs]
			neighbor_q1s, neighbor_q2s = u1u2_to_q1q2(neighbor_lda1s, neighbor_lda2s) ### verified, this works
			neighbor_rp_meters = rp_meters[neighbor_idxs]
			neighbor_a_over_rstars = a_over_rstars[neighbor_idxs]



		### FOR TESTING PURPOSES!
		#if (mplc.rprstar)**2 < 0.05:
		#	continue

		koi_lcfilename = mpdir+'/saved_lcs/'+koi+'_kepler_lightcurve.tsv'
		try:
			koi_lcfile = open(koi_lcfilename, mode='r')
			print('loaded koi_lcfile.')
		except:
			try:
				print('loading lc through moonpy.')
				mplc = MoonpyLC(targetID=koi, clobber='y')
				mplc.get_neighbors(clobber_lc='y') ### produces a light curve file, with all transits, for every planet, ALREADY FLAGGED!
				koi_lcfile = open(koi_lcfilename, mode='r')
			except:
				print('something went wrong getting this light curve through MoonPy.')
				print(' ')
				print(' ')
				continue
		BKJD, fluxes, errors, flags, quarter, in_transit, transiter, target_in_transit, neighbor_in_transit = [], [], [], [], [], [], [], [], []


		for nline,line in enumerate(koi_lcfile):
			linesplit = line.split()
			if nline == 0:
				if len(linesplit) == 7:
					already_detrended = 'n'
				elif len(linesplit) == 9:
					already_detrended = 'y'
				else:
					#raise Exception('something weird with this file.')
					mplc = MoonpyLC(targetID=koi, clobber='y')
					mplc.get_neighbors(clobber_lc='y') ### produces a light curve file, with all transits, for every planet, ALREADY FLAGGED!		
			koi_lcfile.close()
			break

		koi_lcfile = open(koi_lcfilename, mode='r')
		for nline, line in enumerate(koi_lcfile):
			linesplit = line.split()
			if nline == 0:
				if len(linesplit) == 7:
					already_detrended = 'n'
				elif len(linesplit) == 9:
					already_detrended = 'y'

			else:
				linesplit = line.split()
				BKJD.append(float(linesplit[0]))
				fluxes.append(float(linesplit[1]))
				errors.append(float(linesplit[2]))

				if already_detrended == 'n': ### detrending has not occurred!
					flags.append(float(linesplit[3]))
					quarter.append(float(linesplit[4]))
					in_transit.append(linesplit[5])
					try:
						transiter.append(linesplit[6])
					except:
						transiter.append('none')
				
				elif already_detrended == 'y':
					flags.append(float(linesplit[5]))
					quarter.append(float(linesplit[6]))
					in_transit.append(linesplit[7])
					try:
						transiter.append(linesplit[8])
					except:
						transiter.append('none')	

				if in_transit == 'y':
					### identify whether the target, or neighbors are in transit
					all_transiters = line[line.find('[')]
					if koi in all_transiters:
						target_in_transit.append('y')
						neighbor_in_transit.append('n')
					else:
						target_in_transit.append('n')
						neighbor_in_transit.append('y')
				elif in_transit == 'n':
					target_in_transit.append('n')
					neighbor_in_transit.append('n')


		BKJD, fluxes, errors, flags, quarter, in_transit, transiter, target_in_transit, neighbor_in_transit = np.array(BKJD), np.array(fluxes), np.array(errors), np.array(flags), np.array(quarter), np.array(in_transit), np.array(transiter), np.array(target_in_transit), np.array(neighbor_in_transit)

		print('len(BKJD) = ', len(BKJD))
		print('len(quarter) = ', len(quarter))

		unique_quarters = np.unique(quarter)

		### need to expand the transit flagging

		generous_in_transit = []
		for it in in_transit:
			generous_in_transit.append(it)
		generous_in_transit = np.array(generous_in_transit)

		for nit, it in enumerate(in_transit):
			if it == 'y':
				generous_in_transit[nit-5:nit+6] = 'y'

		in_transit = generous_in_transit

		original_in_transit = in_transit #### need this to track each iteration!
		original_transiter = transiter 
		original_target_in_transit = target_in_transit
		original_neighbor_in_transit = neighbor_in_transit 
		
		### CREATE A BATMAN MODEL
		target_batman_params = batman.TransitParams()
		target_batman_params.t0 = koi_tau0
		target_batman_params.per = koi_period ### in days
		target_batman_params.rp = koi_rprstar ### natively in stellar units
		target_batman_params.a = koi_a_over_rstar
		target_batman_params.inc = koi_inc
		target_batman_params.ecc = koi_eccen
		target_batman_params.w = koi_longp
		#u1, u2 = q1q2_to_u1u2(q1, q2)
		target_batman_params.u = [koi_lda1, koi_lda2]
		target_batman_params.limb_dark = 'quadratic'

		target_batman_model = batman.TransitModel(target_batman_params, BKJD)
		target_batman_fluxes = target_batman_model.light_curve(target_batman_params)

		target_delta_fluxes = 1 - target_batman_fluxes

		if len(neighbor_kepois) != 0:
			### make a flux model for every neighbor!
			delta_fluxes = np.linspace(0,0,len(BKJD)) ### start out with zeros!
			for i in np.arange(0,len(neighbor_kepois),1):
				neighbor_batman_params = batman.TransitParams()
				neighbor_batman_params.t0 = neighbor_tau0s[i]
				neighbor_batman_params.per = neighbor_periods[i] ### in days
				neighbor_batman_params.rp = neighbor_rprstars[i] ### natively in stellar units
				neighbor_batman_params.a = neighbor_a_over_rstars[i]
				neighbor_batman_params.inc = neighbor_incs[i]
				neighbor_batman_params.ecc = neighbor_eccens[i]
				neighbor_batman_params.w = neighbor_longps[i]
				#u1, u2 = q1q2_to_u1u2(q1, q2)
				neighbor_batman_params.u = [neighbor_lda1s[i], neighbor_lda2s[i]]
				neighbor_batman_params.limb_dark = 'quadratic'

				neighbor_batman_model = batman.TransitModel(neighbor_batman_params, BKJD)
				neighbor_batman_fluxes = neighbor_batman_model.light_curve(neighbor_batman_params)	
				neighbor_delta_fluxes = 1 - neighbor_batman_fluxes ### out-of-transit = 0, in-transit will be some value (say, 0.1)
				delta_fluxes = delta_fluxes + neighbor_delta_fluxes ### these are what will get subtracted from 1!


			all_neighbor_fluxes = 1 - delta_fluxes 
			target_and_neighbor_fluxes = 1 - delta_fluxes - target_delta_fluxes



		### NOW, RUN A MEDIAN FILTER THROUGH THE LIGHT CURVE
		### WHATEVER IS MARKED IN TRANSIT, MASK IT, AND GENERATE REALISTIC NOISE
		nonmasked_idxs = np.where(in_transit == 'n')[0]
		masked_idxs = np.where(in_transit == 'y')[0]

		target_nonmasked_idxs = np.where(target_in_transit == 'n')[0]
		target_masked_idxs = np.where(target_in_transit == 'y')[0]

		neighbor_nonmasked_idxs = np.where(neighbor_in_transit == 'n')[0]
		neighbor_masked_idxs = np.where(neighbor_in_transit == 'y')[0]

		print('len(BKJD) = ', len(BKJD))
		print('len(nonmasked_idxs) = ', len(nonmasked_idxs))
		print('len(masked_idxs) = ', len(masked_idxs))

		"""
		#NOVEMBER 8th, 2019 -- GPs are way too slow, fall over too much, and don't fix the problem of having locally flat
		#trend lines through masked transits.

		try:
			print('fitting a GP.')
			### DO IT QUARTER BY QUARTER
			final_gp_fluxes, final_gp_errors = [], []

			for quart in unique_quarters:
				print("quarter = ", quart)
				quarter_idxs = np.where(quarter == quart)[0]
				quarter_nonmasked = []
				for nqi, qi in enumerate(quarter_idxs):
					if qi in nonmasked_idxs:
						quarter_nonmasked.append(nqi) ### need to be counting from zero!
				quarter_nonmasked_idxs = np.array(quarter_nonmasked)

				quarter_BKJD, quarter_fluxes, quarter_errors = BKJD[quarter_idxs], fluxes[quarter_idxs], errors[quarter_idxs]
				quarter_gp_fluxes, quarter_gp_errors = GPfit(quarter_BKJD, quarter_fluxes, quarter_errors, quarter_nonmasked_idxs)
				final_gp_fluxes, final_gp_errors = np.concatenate((final_gp_fluxes, quarter_gp_fluxes)), np.concatenate((final_gp_errors, quarter_gp_errors))

			lc_noise = np.nanstd(fluxes[nonmasked_idxs] - final_gp_fluxes[nonmasked_idxs])

		"""


		#except:
		#try:
		print('using a median filter.')
		nonmasked_medfilt = medfilt(fluxes[nonmasked_idxs], kernel_size=9)
		### COMPUTE THE STANDARD DEVIATION OF THE NON-MASKED FLUXES FROM THE MEDIAN FILTER
		lc_noise = np.nanstd(fluxes[nonmasked_idxs] - nonmasked_medfilt)

		medfilt_interpolator = interp1d(BKJD[nonmasked_idxs], nonmasked_medfilt, bounds_error=False, fill_value='extrapolate')
		interpolated_median_filter = medfilt_interpolator(BKJD) ### this is the median filter at all time stamps!	
		#except:
		#	traceback.print_exc()
		#	continue 


		trend_line = interpolated_median_filter


		#### USING A MEDIAN FILTER WHETHER WE RANDOMIZE THE TRANSIT LOCATION OR NOT!
		### November 8th, 2019 -- the median filter is the cleanest way to clean all transits.
		transit_removed_fluxes = []
		for flux in fluxes:
			transit_removed_fluxes.append(flux) ### THIS IS JUST SO TRANSIT_REMOVED_FLUX IS NOT POINTING TO FLUXES.

		transit_removed_fluxes = np.array(transit_removed_fluxes)
		transit_removed_fluxes[masked_idxs] = np.random.normal(loc=trend_line[masked_idxs], scale=lc_noise)

		### how many transits were observed during the actual Kepler mission
		ntransits_per_one_mission = (np.nanmax(BKJD) - np.nanmin(BKJD)) / koi_period

		### how many times do we need to iterate on this?
		n_injection_iterations = ntransits_per_star / ntransits_per_one_mission ### if we want 100, and we saw max 5, we need 20 injection iterations!

		if n_injection_iterations < 1:
			n_injection_iterations = 1
		elif n_injection_iterations > ntransits_per_star: ### this shouldn't be the case!
			continue 


		for ii in np.arange(0,n_injection_iterations,1):
			### for each iteration of injection!
			iteration_number = ii 
			print('injection iteration = ', iteration_number)

			fp_filename = koi+'_kepler_fptest_lightcurve_iteration'+str(int(iteration_number))+'.tsv'
			if randomize_injection_location == 'n':
				fp_file_destination = '/Users/hal9000/Documents/Software/MoonPy/fptest_planet_injections'
				
			elif randomize_injection_location == 'y':
				fp_file_destination = '/Users/hal9000/Documents/Software/MoonPy/fptest_planet_injections_random_loc'


			### check if this file already exists! if it does, skip it!
			if (os.path.exists(fp_file_destination+'/'+fp_filename)) and (overwrite_fplc == 'n'):
				print('file already exists! SKIPPING.')
				break



			if randomize_injection_location == 'n':
				### IN THIS CASE, WE'RE GOING TO MASK OUT ALL TRANSITS, USE A MEDIAN FILTER, AND INTERPOLATE THROUGH THE BASELINE.
				### place the injected transits right where the planets are right now.
				injected_transit_fluxes = transit_removed_fluxes*target_batman_fluxes
				new_tau0 = koi_tau0

			elif randomize_injection_location == 'y':

				### IN THIS CASE, WE WANT TO DIVIDE OUT THE BATMAN MODEL, LEAVING RESIDUAL CORRELATED NOISE!
				#transit_removed_fluxes = fluxes / target_batman_fluxes ### DIVIDING OUT THE MODEL! Leaves in tact the correlated noise.

				### for drawing a line where the transits are.
				#target_batman_model_times_masked_filter = trend_line*target_batman_fluxes

				### divide out the neighbors!
				#if len(neighbor_kepois) != 0:
				#	all_transits_model_times_masked_filter = trend_line*target_and_neighbor_fluxes
				#	#### THIS REMOVES ALL THE NEIGHBOR TRANSITS -- further modifies transit_removed_fluxes from above.
				#	transit_removed_fluxes = transit_removed_fluxes / all_neighbor_fluxes 

				### we want to produce 100 transits, regardless of planet period!
				### do this WITHOUT MANIPULATING THE PLANET'S PERIOD (THAT AFFECTS THE TRANSIT DURATION).

				### choose some random index ### PLACE 
				random_start_idx = np.random.randint(low=0, high=len(BKJD))
				tau0_offset = BKJD[random_start_idx] - BKJD[0] ### this is how much you need to shift tau0!!!!
				### ^^^ IMPORTANT FOR HAVING THE RIGHT EPHEMERIS IN THE CNNLC PARSER!
				new_tau0 = koi_tau0 + tau0_offset

				shifted_batman_fluxes = np.concatenate((target_batman_fluxes[random_start_idx:], target_batman_fluxes[:random_start_idx]))
				injected_transit_fluxes = transit_removed_fluxes*shifted_batman_fluxes
				in_transit = np.concatenate((original_in_transit[random_start_idx:], original_in_transit[:random_start_idx])) ### we're shifting around the in_transit locations, because we've shifted the transit locations!
				transiter = np.concatenate((original_transiter[random_start_idx:], original_transiter[:random_start_idx])) ### ditto above.

				new_in_transit_idxs = np.where(in_transit == 'y')[0]

				### NOTE THAT in_transit will be used for parsing the CNNLCs (fptest_cnnlc_generator) for cutting up these shifted light curves!

			### test that this is working!
			if show_plots == 'y':
				fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
				ax1.scatter(BKJD, fluxes, c='LightCoral', s=10)
				#if randomize_injection_location == 'n':
				ax1.plot(BKJD, trend_line, c='k', linewidth=1)
				#elif randomize_injection_location == 'y':
					#if len(neighbor_kepois) == 0:
					#	ax1.plot(BKJD, target_batman_model_times_masked_filter, c='k', linewidth=1)
					#else:
					#	ax1.plot(BKJD, all_transits_model_times_masked_filter, c='k', linewidth=1)

				ax1.set_ylabel('Real')

				ax1.scatter(BKJD[masked_idxs], fluxes[masked_idxs], c='g', s=10)
				ax2.scatter(BKJD, transit_removed_fluxes, c='DodgerBlue', s=10)
				if randomize_injection_location == 'n':
					ax2.plot(BKJD, trend_line, c='k', linewidth=1)
				ax2.set_ylabel('transits removed')

				ax3.scatter(BKJD, injected_transit_fluxes, c='k', s=10)
				ax3.scatter(BKJD[new_in_transit_idxs], injected_transit_fluxes[new_in_transit_idxs], c='g', s=10)

				ax3.set_ylabel('transits injected')

				ax1.set_title(koi+r', $\sigma$ = '+str(round(lc_noise, 2)))
				plt.show()


			if write_to_file == 'y':
				fp_file = open(fp_file_destination+'/'+fp_filename, mode='w')	

				fp_file.write('new_tau0 = '+str(new_tau0)+'\n')
				fp_file.write('BKJD\tfluxes\terrors\tflags\tquarter\tin_transit\ttransiter\n')
				for bkjd, flux, error, flag, quart, tl, pl in zip(BKJD, injected_transit_fluxes, errors, flags, quarter, in_transit, transiter):
					fp_file.write(str(bkjd)+'\t'+str(flux)+'\t'+str(error)+'\t'+str(quart)+'\t'+str(tl)+'\t'+str(pl)+'\n')
				fp_file.close()
				print('false positive test file written.')

			print('light curve fully processed.')
			print(" ")
			print(" ")

	except:
		traceback.print_exc()
		time.sleep(5)
		raise Exception('triage the problem above.')
		print 




